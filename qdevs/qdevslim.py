
"""Generic LIM and QDEVS-LIM Model."""

from qdevs import *

_INF = float("inf")
_EPS = 1e-15

class QdevsLimNode(QdevsDevice):

    def __init__(self, C, R, I, granularity=None, v0=0.0):

        """Implements a generic LIM node ODE.
        """

        QdevsDevice.__init__(self, granularity, v0)

        self.R = R
        self.C = C
        self.I = I
        self.branch_connections = []

    def connect_branch(self, branch, polarity):

        self.branch_connections.append((branch, polarity))
        
        self.connect_outputs(branch)

    def evaluate(self, time):

        self.tnext = _INF

        isum = 0.0

        for branch, polarity in self.branch_connections:
            isum += polarity * branch.state

        #denom = self.I * self.R - self.state - self.R * isum
        denom = 2.0 * self.I * self.R - 2.0 * self.state - 2.0 * self.R * isum - self.granularity


        if denom != 0.0:

            #delta_time = self.C * self.R * self.granularity / denom
            delta_time = 2.0 * self.C * self.R * self.granularity / denom

            if delta_time > 0.0:
                self.tnext = self.tlast + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = self.tlast - delta_time
                self.trajectory_direction = -1.0

            else:
                self.tnext = self.tlast + _EPS


class QdevsLimBranch(QdevsDevice):

    """Implements a two-port QDEVS LIM Branch ODE.
    """

    def __init__(self, L, R, V, granularity=None, i0=0.0):

        QdevsDevice.__init__(self, granularity, i0)

        self.R = R
        self.L = L
        self.V = V

        self.nodei = None
        self.nodej = None

    def connect_nodes(self, nodei, nodej):

        self.nodei = nodei
        self.nodej = nodej

        self.connect_outputs(nodei, nodej)

        nodei.connect_branch(self, polarity=1.0)
        nodej.connect_branch(self, polarity=-1.0)

    def evaluate(self, time):

        self.tnext = _INF

        vij = self.nodei.state - self.nodej.state

        #denom = self.V + vij - self.R * self.state
        denom = 2.0 * vij - 2.0 * self.R * self.state - self.R * self.granularity + 2.0 * self.V

        if denom != 0.0:

            #delta_time = self.L * self.granularity / denom
            delta_time = 2.0 * self.L * self.granularity / denom

            if delta_time > 0.0:
                self.tnext = self.tlast + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = self.tlast + -delta_time
                self.trajectory_direction = -1.0

            else:
                self.tnext = self.tlast + _EPS


class QdevsLimGround(QdevsDevice):

    """ LIM ground. LimBranchs can be connected to this.
    """

    def __init__(self):

        QdevsDevice.__init__(self)

    def connect_branch(self, branch, polarity):

        pass

    def evaluate(self, time):

        pass


class QdevsLimSystem(DevsSystem):

    def __init__(self, current_granularity=1e-3, voltage_granularity=1e-3):
        
        DevsSystem.__init__(self)
        self.current_granularity = current_granularity
        self.voltage_granularity = voltage_granularity
        self.nodes = []
        self.branches = []

    def add_node(self, C, R=None, I=None, granularity=None):

        """Adds a node shunt to the system
        """

        device = QdevsLimNode(C, R, I)
        if granularity:
            device.granularity = granularity
        else:
            device.granularity = self.voltage_granularity
        self.nodes.append(device)
        self.devices.append(device)

        return device

    def add_branch(self, nodei, nodej, L, R=None, V=None, granularity=None):

        """Adds a branch to the system
        """

        device = QdevsLimBranch(L, R, V)
        if granularity:
            device.granularity = granularity
        else:
            device.granularity = self.current_granularity
        self.branches.append(device)
        self.devices.append(device)
        device.connect_nodes(nodei, nodej)

        return device

    def add_ground(self):

        device = QdevsLimGround()
        self.nodes.append(device)
        self.devices.append(device)

        return device

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

