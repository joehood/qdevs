
"""Generic QDEVS  models and simulator with support for LIM
branch/node electrical circuit representation."""


import numpy as np
import numpy.linalg as la

from qdevs import *

_INF = float("inf")
_EPS = 1e-15


class StateSpace(object):

    def __init__(self, a, b, c=None, d=None, x0=None, u0=None):

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.x0 = x0
        self.x = x0
        self.u0 = u0
        self.n = self.a.shape[0]
        self.m = self.b.shape[1]
        self.dt = -1.0

    def initialize(self, dt, u0=None):

        self.dt = dt

        self.n = self.a.shape[0]
        self.m = self.b.shape[1]

        if self.c is None:
            self.c = np.eye(self.n)

        self.p = self.c.shape[0]

        if self.d is None:
            self.d = np.zeros((self.p, self.m))

        if self.x0 is None:
            self.x = np.zeros((self.n, 1))
            
        else:
            self.x = self.x0

        if u0 is not None:
            self.u = u0
        elif self.u0 is not None:
            self.u = self.u0
        else:
            self.u = np.zeros((self.m, 1))

        eye = np.eye(self.n)

        self.apr = la.inv(eye - dt * self.a)
        self.bpr = np.dot(self.apr, dt * self.b)

        self.y = np.dot(self.c, self.x) + np.dot(self.d, self.u)

        return self.y

    def step(self, u):

        self.u = u
        self.x = np.dot(self.apr, self.x) + np.dot(self.bpr, self.u)
        self.y = np.dot(self.c, self.x) + np.dot(self.d, self.u)

        return self.y

    def run(self, tf):

        t = np.arange(0.0, tf, self.dt)
        y = np.zeros((self.n, t.size))
        y[:,0:1] = self.y

        for i in range(1, t.size):
            y[:,i:i+1] = self.step(self.u0)

        return t, y


class QdevsLimNode(QdevsDevice):

    """Implements a generic LIM node shunt ODE.
    """

    def __init__(self, C, R, I, v0=0.0, granularity=None):

        QdevsDevice.__init__(self, v0, granularity)

        self.R = R
        self.C = C
        self.I = I

        self.index = -1
        self.branch_connections = []

    def connect_branch(self, branch, polarity):

        """Adds a branch to the connected branch list with the
        given polarity, and adds the branch as an output device
        so external transitions will be triggered at the branch
        when this node is updated.
        """

        self.branch_connections.append((branch, polarity))
        
        self.connect_outputs(branch)

    def update(self, time):

        """Calculates the next time to reach the next quantized
        state using a semi-implicit integration.
        """

        isum = 0.0

        for branch, polarity in self.branch_connections:
            isum += polarity * branch.state

        self.last_state = self.state
        dt = time - self.tlast
        next_dt = _INF

        self.internal_state += self.derivative * dt

        if self.internal_state >= self.state + self.granularity - self.epsilon:
            self.state += self.granularity

        elif self.internal_state <= self.state - 0.5 * self.granularity + self.epsilon:
            self.state -= self.granularity

        self.derivative = -(self.state + self.R * isum - self.I * self.R) / (self.C * self.R)

        if self.derivative > 0.0:
            next_dt = (self.state + self.granularity - self.internal_state) / self.derivative
        
        elif self.derivative < 0.0:
            next_dt = (self.state - 0.5 * self.granularity - self.internal_state) / self.derivative 
        
        self.tnext = time + next_dt
        self.tlast = time
        self.save(time)


class QdevsLimBranch(QdevsDevice):

    """Implements a two-port QDEVS LIM Branch ODE.
    """

    def __init__(self, L, R, V, i0=0.0, granularity=None):

        QdevsDevice.__init__(self, i0, granularity)

        self.R = R
        self.L = L
        self.V = V

        self.index = -1
        self.nodei = None
        self.nodej = None

    def connect_nodes(self, nodei, nodej):

        """Specifies the ith and jth node devices to which this
        branch is connected so we can get the updated node voltages
        needed for evaluate. Also adds these nodes as output devices
        so external transitions will be triggered at the nodes when this
        branch's state is updated.
        """

        self.nodei = nodei
        self.nodej = nodej

        self.connect_outputs(nodei, nodej)

        nodei.connect_branch(self, polarity=1.0)
        nodej.connect_branch(self, polarity=-1.0)

    def update(self, time):

        """Calculates the next time to reach the next quantized
        state using a semi-implicit integration.
        """

        vij = self.nodei.state - self.nodej.state

        self.last_state = self.state
        dt = time - self.tlast
        next_dt = _INF

        self.internal_state += self.derivative * dt

        if self.internal_state >= self.state + self.granularity - self.epsilon:
            self.state += self.granularity

        elif self.internal_state <= self.state - 0.5 * self.granularity + self.epsilon:
            self.state -= self.granularity

        self.derivative = vij - self.R * self.state + self.V / self.L

        if self.derivative > 0.0:
            next_dt = (self.state + self.granularity - self.internal_state) / self.derivative
        
        elif self.derivative < 0.0:
            next_dt = (self.state - 0.5 * self.granularity - self.internal_state) / self.derivative 
        
        self.tnext = time + next_dt
        self.tlast = time
        self.save(time)


class QdevsLimGround(QdevsDevice):

    """ LIM ground. LimBranchs can be connected to this. It
    forces a constant reference voltage.
    """

    def __init__(self, vref=0.0):

        QdevsDevice.__init__(self)

        self.state0 = vref
        self.state = vref
        self.index = 0

    def connect_branch(self, branch, polarity):

        pass

    def update(self, time):

        pass


class QdevsLimSystem(DevsSystem):

    """A QDEVS lim system representation and simulator that with some
    additional automations for LIM-based devices and circuit topologies.
    """

    def __init__(self, current_granularity=1e-3, voltage_granularity=1e-3):
        
        DevsSystem.__init__(self)
        self.current_granularity = current_granularity
        self.voltage_granularity = voltage_granularity
        self.nodes = []
        self.branches = []
        self.ss = None
        self.node_index = 1
        self.branch_index = 0

    def build_ss(self):

        n = len(self.nodes)
        m = len(self.branches)

        a = np.zeros((n+m, n+m))
        b = np.zeros((n+m, n+m))
        u = np.zeros((n+m, 1))

        for node in self.nodes:

            ii = node.index

            if ii == 0:
                continue

            a[ii, ii] = -1.0 / (node.R * node.C)
            b[ii, ii] = 1.0 / node.C
            u[ii, 0] = node.I

            for branch in self.branches:

                k = branch.index
                i = branch.nodei.index
                j = branch.nodej.index

                if ii == i:
                    a[i, n+k] = -1.0 / node.C

                elif ii == j:
                    a[j, n+k] = 1.0 / node.C

        for branch in self.branches:

            k = branch.index
            i = branch.nodei.index
            j = branch.nodej.index

            a[n+k, n+k] = -branch.R / branch.L
            a[n+k, i] = 1.0 / branch.L
            a[n+k, j] = -1.0 / branch.L

            b[n+k, n+k] = 1.0 / branch.L
            u[n+k, 0] = branch.V

        self.ss = StateSpace(a[1:, 1:], b[1:, 1:], u0=u[1:,:])

    def initialize(self, t0=0.0):

        self.build_ss()

        super(QdevsLimSystem, self).initialize(t0)

    def add_node(self, C, R=None, I=None, granularity=None):

        """Adds a node shunt to the system.
        """

        device = QdevsLimNode(C, R, I)
        device.index = self.node_index
        if granularity:
            device.granularity = granularity
        else:
            device.granularity = self.voltage_granularity
        self.nodes.append(device)
        self.devices.append(device)

        self.node_index += 1

        return device

    def add_branch(self, nodei, nodej, L, R=None, V=None, granularity=None):

        """Adds a branch to the system.
        """

        device = QdevsLimBranch(L, R, V)

        device.index = self.branch_index

        if granularity:
            device.granularity = granularity
        else:
            device.granularity = self.current_granularity

        self.branches.append(device)
        self.devices.append(device)

        device.connect_nodes(nodei, nodej)

        self.branch_index += 1

        return device

    def add_ground(self, vref=0.0):

        """Adds a ground (reference node) to the system.
        """

        device = QdevsLimGround(vref)
        self.nodes.append(device)
        self.devices.append(device)

        return device

    def advance(self):

        """Advances the simulation to the next scheduled event.
        Currently this is the same as the QdevsSystem advance routine,
        but new behavior and optimizations can be added here for
        exploiting features, stabilizing and control errors specific
        to LIM-based topologies.
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


