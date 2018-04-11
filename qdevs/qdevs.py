"""Generic QDEVS Model."""


class DevsDevice(object):

    """Atomic QDEVS Device"""

    def __init__(self, state0=0.0):

        self.state = 0.0
        self.state0 = 0.0

        self.tnext = float("inf")
        self.tlast = 0.0

        self.input = 0.0

        self.output_devices = []
        self.input_devices = []

        self.time_history = []
        self.state_history = []

    def connect_outputs(self, *devices):

        for device in devices:
            self.output_devices.append(device)
            device.input_devices.append(self)

    def broadcast(self, time):

        for output_device in self.output_devices:
            output_device.external_transition(self, time, self.state)

    def save(self, time, reset=False):

        if reset:
            self.time_history = [time]
            self.state_history = [self.state]
        else:
            self.time_history.append(time)
            self.state_history.append(self.state)

    def external_transition(self, sender, time, value):

        self.input = value
        self.evaluate(time)

    def initialize(self, time):

        throw(NotImplementedError)

    def internal_transition(self, time):

        throw(NotImplementedError)

    def evaluate(self, time):

        throw(NotImplementedError)


class QdevsDevice(DevsDevice):

    """Atomic QDEVS Device"""

    def __init__(self, granularity=1e-3, state0=0.0):

        DevsDevice.__init__(self, state0)

        self.granularity = granularity
        self.trajectory_direction = 0.0

    def initialize(self, time):

        self.state = self.state0
        self.tnext = float("inf")
        self.tlast = time
        self.save(time)
        self.evaluate(time)

    def internal_transition(self, time):

        self.state += self.trajectory_direction * self.granularity
        self.tlast = time
        self.save(time)
        self.evaluate(time)


class QdevsSystem(object):

    def __init__(self):

        self.devices = []
        self.time = 0.0

    def add_devices(self, *devices):

        for device in devices:
            self.devices.append(device)

    def initialize(self, t0=0.0):

        self.time = t0

        for device in self.devices:
            device.initialize(t0)

        for device in self.devices:
            device.broadcast(t0)

    def run(self, tstop):

        while(self.time < tstop):
            self.advance()

    def advance(self):

        tnext = float("inf")

        for device in self.devices:
            tnext = min(tnext, device.tnext)

        self.time = tnext

        imminent_devices = []

        for device in self.devices:
            if device.tnext <= self.time:
                imminent_devices.append(device)

        for device in imminent_devices:
            device.internal_transition(self.time)

        for device in imminent_devices:
            device.broadcast(self.time)




