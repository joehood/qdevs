"""Generic DEVS and QDEVS Models."""

from __future__ import division

from collections import deque

from scipy.signal import resample
from scipy.interpolate import interp1d
import numpy as np

_INF = float("inf")
_EPS = 1e-9

class DevsEvent(object):

    """Generic DEVS Event"""

    def __init__(self, sender, time, value):

        self.sender = sender
        self.time = time
        self.value = value


class DevsDevice(object):

    """Generic Atomic DEVS Device."""

    def __init__(self, state0=0.0):

        self.state0 = state0
        self.state = state0
        self.last_state = state0
        self.tnext = _INF
        self.tlast = 0.0
        self.input = 0.0
        self.sender = None
        self.input_events = deque()
        self.output_devices = []
        self.time_history = []
        self.state_history = []

    def connect_outputs(self, *devices):

        """Connect this device to an output devices. When this
        device goes through an internal transistion, it will trigger
        an external event on these devices and send the event data.
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
            event = self.input_events.pop()
            self.sender = event.sender
            self.input = event.value
            self.update(event.time)

    def broadcast(self, time):

        """Send external events to the connected output devices.
        """

        if self.state != self.last_state:
            for output_device in self.output_devices:
                output_device.add_input(DevsEvent(self, time, self.state))

    def save(self, time, reset=False):

        """Save the current time and state to the history arrays.
        """

        if reset:
            self.time_history = [time]
            self.state_history = [self.state]

        elif self.state != self.last_state:
            self.time_history.append(time)
            self.state_history.append(self.state)

    def initialize(self, time):

        """Can be overridden in derived class. This is called at the
        beginning of the simulation. Usually, initial states and the
        initial tnext values are set here.
        """

        self.state = self.state0
        self.last_state = self.state0
        self.tlast = time
        self.tnext = _INF
        self.save(time, reset=True)
        self.broadcast(time)

    def update(self, time):

        """Must be implemented in derived class. This will be called when
        the simulation advances to the current tnext value of this
        device. Usually the state is updated to the appropriate next
        value here.
        """

        raise NotImplementedError()


class QdevsDevice(DevsDevice):

    """Generic Atomic QDEVS Device. Contains some additional data and
     functionality speicific to Quantized DEVS devices.
    """

    def __init__(self, state0=0.0, granularity=1e-3, epsilon=None):

        DevsDevice.__init__(self, state0)

        self.granularity = granularity

        if epsilon:
            self.epsilon = epsilon
        elif granularity:
            self.epsilon = 0.5 * granularity

        self.internal_state = state0
        self.derivative = 0.0
        self.epsilon = 0.0
        
    def initialize(self, time):

        self.state = self.state0
        self.internal_state = self.state0
        self.derivative = 0.0
        self.tlast = time
        self.tnext = _INF
        self.update(time)
        self.save(time, reset=True)
        self.broadcast(time)


class DevsSystem(object):

    """Generic DEVS system representation and simulator."""

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
        initialize() must be called before running for the first time.
        The simulator can have multiple run() calls in the same
        simulation to enable external events to be implemented.
        """

        self.tstop = tstop

        while(self.time < tstop):
            self.advance()

    def advance(self):

        """Advances the simulation to the next scheduled event, 
        imminent devices will have internal transitions, those devices
        will broadcasts events to their output devices who will then
        process those events.
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
            device.update(self.time)

        for device in imminent_devices:
            device.broadcast(self.time)

        for device in self.devices:
            device.process_inputs()


class QdevsSystem(DevsSystem):

    """ Generic QDEVS system representation and simulator. Contains
    specific additions for handling quantized devices.
    """

    def __init__(self, granularity=1e-3, epsilon=None):
        
        DevsSystem.__init__(self)

        if granularity:
            self.granularity = granularity
        else:
            self.granularity =  1e-3

        if epsilon:
            self.epsilon = epsilon
        else:
            self.epsilon = 0.25 * granularity

    def add_devices(self, *devices):

        """Adds one or more devices to the system and cascades
        the default granularity to the devices.
        """

        for device in devices:

            if isinstance(device, QdevsDevice):

                if not device.granularity:
                    device.granularity = self.granularity

                if not device.epsilon:
                    device.epsilon = self.epsilon

            self.devices.append(device)


class ConstantSource(DevsDevice):

    """Constant source model
    """

    def __init__(self, value):

        DevsDevice.__init__(self)

        self.state0 = value
        self.state = value

    def initialize(self, time):

        self.state = self.state0
        self.tlast = time
        self.tnext = _INF
        self.save(time, reset=True)
        self.broadcast(time)

    def set_value(self, value):

        if value != self.value:
            self.value = value
        self.tnext = self.tlast

    def update(self, time):

        self.save(time)


class SquareWaveSource(DevsDevice):

    """Simple square wave with variable duty and zero rise/fall
    time.
    """

    def __init__(self, x1, x2, t1, t2):

        DevsDevice.__init__(self)

        self.state0 = x1
        self.x1 = x1
        self.x2 = x2
        self.t1 = t1
        self.t2 = t2

    def initialize(self, time):

        self.state = self.state0
        self.tlast = time
        self.tnext = time + self.t1
        self.save(time, reset=True)
        self.broadcast(time)

    def update(self, time):

        self.last_state = self.state
        self.tnext = self.t1 + self.t2

        if self.state == self.x1:
            self.state = self.x2
            self.tnext = time + self.t2

        elif self.state == self.x2:
            self.state = self.x1
            self.tnext = time + self.t1

        self.tlast = time
        self.save(time)
        

class Integrator(QdevsDevice):

    """Simple linear integrator with gain and no limits with form:
    x' = k*u
    """

    def __init__(self, gain, x0=0.0, granularity=None, epsilon=None):

        QdevsDevice.__init__(self, x0, granularity, epsilon)

        self.gain = gain

    def update(self, time):

        self.last_state = self.state
        dt = time - self.tlast
        next_dt = _INF

        self.internal_state += self.derivative * dt

        if self.internal_state >= self.state + self.granularity - self.epsilon:
            self.state += self.granularity
            self.broadcast(time)

        elif self.internal_state <= self.state - 0.5 * self.granularity + self.epsilon:
            self.state -= self.granularity
            self.broadcast(time)

        self.derivative = self.gain * self.input

        if self.derivative > 0.0:
            next_dt = (self.state + self.granularity - self.internal_state) / self.derivative
        
        elif self.derivative < 0.0:
            next_dt = (self.state - 0.5 * self.granularity - self.internal_state) / self.derivative 
        
        self.tnext = time + abs(next_dt)
        self.tlast = time
        self.save(time)
        

class DifferentialEquation(QdevsDevice):

    """Represents a continuous first order ODE of the form:
    x' = a * x + b * u
    """

    def __init__(self, a, b, x0=0.0, granularity=None, epsilon=None):

        QdevsDevice.__init__(self, x0, granularity, epsilon)

        self.a = a
        self.b = b

    def update(self, time):

        self.last_state = self.state
        dt = time - self.tlast
        next_dt = _INF

        self.internal_state += self.derivative * dt

        if self.internal_state >= self.state + self.granularity - self.epsilon:
            self.state += self.granularity
            self.broadcast(time)

        elif self.internal_state <= self.state - 0.5 * self.granularity + self.epsilon:
            self.state -= self.granularity
            self.broadcast(time)

        self.derivative = self.a * self.internal_state + self.b * self.input

        if self.derivative > 0.0:
            next_dt = (self.state + self.granularity - self.internal_state) / self.derivative
        
        elif self.derivative < 0.0:
            next_dt = (self.state - 0.5 * self.granularity - self.internal_state) / self.derivative 
        
        self.tnext = time + abs(next_dt)
        self.tlast = time
        self.save(time)
        

def resample(times, values, tf, npoints=1000):

    """Resamples the given time/value event arrays from time 0 to tf
    for npoints using a zero-order hold. This is useful for plotting
    results and quantifying error.
    """

    values.append(values[-1])
    times.append(tf)
    f = interp1d(times, values, kind='zero')
    times2 = np.linspace(times[0], times[-1], npoints)
    values2 = f(times2)
    return times2, values2


