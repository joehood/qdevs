"""Generic DEVS and QDEVS Models."""

from __future__ import division
from collections import deque

_SCIPY = True

try:
    from scipy.signal import resample
    from scipy.interpolate import interp1d
    import numpy as np
except:
    _SCIPY = False

_INF = float("inf")
_EPS = 0.0


class DevsEvent(object):

    def __init__(self, sender, time, value):

        self.sender = sender
        self.time = time
        self.value = value


class DevsDevice(object):

    """Generic Atomic DEVS Device."""

    def __init__(self, state0=0.0):

        self.state = 0.0
        self.state0 = 0.0
        self.tnext = _INF
        self.tlast = 0.0
        self.input = 0.0
        self.input_events = deque()
        self.output_devices = []
        self.time_history = []
        self.state_history = []

    def connect_outputs(self, *devices):

        """Connect this device to aoutput devices. When this device
        goes through an internal transistion, it will trigger
        an external event on those devices and send the event.
        """

        for device in devices:
            self.output_devices.append(device)

    def add_input(self, event):

        """Append an event to the input event queue.
        """

        self.input_events.appendleft(event)

    def process_inputs(self):

        """Processes all input events in the input queue.
        """

        while self.input_events:
            self.external_transition(self.input_events.pop())

    def broadcast(self, time):

        """Trigger an external event on the connected output devices
        and send the event to those devices.
        """

        for output_device in self.output_devices:
            output_device.add_input(DevsEvent(self, time, self.state))

    def save(self, time, reset=False):

        """Save the current time and state to the history arrays.
        """

        if reset:
            self.time_history = [time]
            self.state_history = [self.state]
        else:
            self.time_history.append(time)
            self.state_history.append(self.state)

    def external_transition(self, event):

        """An external input device triggers this provides the event
        data.
        """

        self.input = event.value
        self.evaluate(event.time)

    def initialize(self, time):

        """Must be implemented in derived class. This is called at the
        beginning of the simulation. Usually, initial states and tnext
        values are set here.
        """

        raise NotImplementedError()

    def internal_transition(self, time):

        """Must be implemented in derived class. This will be called when
        the simulation advances to the current tnext value of this
        device. Usually the state is updated to the appropriate next
        value here.
        """

        raise NotImplementedError()

    def evaluate(self, time):

        """Must be implemented in derived class. This is called when it
        is necessary for this device to determine it's next transition
        time (tnext)."""

        raise NotImplementedError()


class QdevsDevice(DevsDevice):

    """Generic Atomic QDEVS Device. Contains some additional data and
     functionality speicific to Quantized DEVS devices.
    """

    def __init__(self, granularity=1e-3, state0=0.0):

        DevsDevice.__init__(self, state0)

        self.granularity = granularity
        self.trajectory_direction = 0.0

    def initialize(self, time):

        self.state = self.state0
        self.trajectory_direction = 0.0
        self.internal_transition(time)

    def internal_transition(self, time):

        self.state += self.trajectory_direction * self.granularity
        self.save(time)
        self.evaluate(time)
        self.tlast = time


class DevsSystem(object):

    def __init__(self):

        self.devices = []
        self.time = 0.0

    def add_devices(self, *devices):

        """Adds one or more devices to the system.
        """

        for device in devices:
            self.devices.append(device)

    def initialize(self, t0=0.0):

        """This should be called at the start of a simulation.
        """

        self.time = t0
        self.tstop = 0.0

        for device in self.devices:
            device.initialize(t0)

        for device in self.devices:
            device.broadcast(t0)

        for device in self.devices:
            device.process_inputs()

    def run(self, tstop):

        """Run the simulation from the current time until tstop. 
        Initialize must be called before running for the fisrt time.
        The simulator can have multiple run calls in the same simulation
        to enable external events to be implemented.
        """
        self.tstop = tstop
        while(self.time < tstop):
            self.advance()

    def advance(self):

        """Advances the simulation to the next scheduled event.
        """

        tnext = _INF

        for device in self.devices:
            tnext = min(tnext, device.tnext)

        self.time = max(tnext, self.time + _EPS)
        
        if self.time > self.tstop:
            return

        imminent_devices = []

        for device in self.devices:
            if device.tnext <= self.time:
                imminent_devices.append(device)

        for device in imminent_devices:
            device.internal_transition(self.time)

        for device in imminent_devices:
            device.broadcast(self.time)

        for device in self.devices:
            device.process_inputs()


class QdevsSystem(DevsSystem):

    def __init__(self, granularity=1e-3):
        
        DevsSystem.__init__(self)
        self.granularity = granularity

    def add_devices(self, *devices):

        """Adds one or more devices to the system.
        """

        for device in devices:
            if isinstance(device, QdevsDevice):
                if not device.granularity:
                    device.granularity = self.granularity
            self.devices.append(device)


class SquareWaveSource(DevsDevice):

    def __init__(self, x1, x2, t1, t2):

        DevsDevice.__init__(self)

        self.x1 = x1
        self.x2 = x2
        self.t1 = t1
        self.t2 = t2

    def initialize(self, time):

        self.tlast = time
        self.state = self.x1
        self.tnext = time + self.t1
        self.save(time, reset=True)

    def internal_transition(self, time):

        if self.state == self.x1:
            self.state = self.x2
            self.tnext = time + self.t2

        elif self.state == self.x2:
            self.state = self.x1
            self.tnext = time + self.t1

        self.save(time)
        self.tlast = time


class Integrator(QdevsDevice):

    def __init__(self, gain, granularity=None, x0=0.0):

        QdevsDevice.__init__(self, granularity, x0)

        self.gain = gain

    def evaluate(self, time):

        self.tnext = float("inf")

        if self.input != 0.0:

            delta_time = self.granularity / (self.gain * self.input)

            if delta_time > 0.0:
                self.tnext = time + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = time - delta_time
                self.trajectory_direction = -1.0
            else:
                self.tnext = time + _EPS


class DifferentialEquation(QdevsDevice):

    def __init__(self, a, b, granularity=None, x0=0.0):

        QdevsDevice.__init__(self, granularity, x0)

        self.a = a
        self.b = b

    def evaluate(self, time):

        self.tnext = float("inf")

        denom = self.a * self.state + self.b * self.input

        if denom != 0.0:

            delta_time = self.granularity / denom

            if delta_time > 0.0:
                self.tnext = time + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = time - delta_time
                self.trajectory_direction = -1.0

            else:
                self.tnext = time + _EPS


def resample(times, values, tf, npoints=1000):

    """Resamples the given time/value event arrays from time 0 to tf
    for npoints using a zero-order hold. This is useful for plotting
    results. Scipy is required.
    """

    if _SCIPY:
        values.append(values[-1])
        times.append(tf)
        f = interp1d(times, values, kind='zero')
        times2 = np.linspace(times[0], times[-1], npoints)
        values2 = f(times2)
        return times2, values2


