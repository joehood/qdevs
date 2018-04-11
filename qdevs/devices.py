"""QDEVS Device Definitions"""

import qdevs


class SquareWaveSource(qdevs.DevsDevice):

    def __init__(self, x1, x2, t1, t2):

        qdevs.DevsDevice.__init__(self)

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


class Integrator(qdevs.QdevsDevice):

    def __init__(self, gain, granularity, x0):

        qdevs.QdevsDevice.__init__(self, granularity, x0)

        self.gain = gain

    def evaluate(self, time):

        self.tnext = float("inf")

        if self.input != 0.0:

            delta_time = self.granularity / (self.gain * self.input)

            if delta_time > 0.0:
                self.tnext = time + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = time + -delta_time
                self.trajectory_direction = -1.0
            else:
                self.tnext = time + delta_time


class DifferetialEquation(qdevs.QdevsDevice):

    def __init__(self, a2, a1, a0, granularity=1e-3, x0=0.0):

        qdevs.QdevsDevice.__init__(self, granularity, x0)

        self.a2 = a2
        self.a1 = a1
        self.a0 = a0

    def evaluate(self, time):

        self.tnext = float("inf")

        denom = self.input + self.a0 - self.a1 * self.state

        if denom != 0.0:

            delta_time = self.a2 * self.granularity / denom

            if delta_time > 0.0:
                self.tnext = time + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = time + -delta_time
                self.trajectory_direction = -1.0

            else:
                self.tnext = time + 1e-15


class LimNode(qdevs.QdevsDevice):

    def __init__(self, C, R, I, granularity=1e-3, v0=0.0):

        qdevs.QdevsDevice.__init__(self, granularity, v0)

        self.R = R
        self.C = C
        self.I = I

    def evaluate(self, time):

        self.tnext = float("inf")

        isum = 0.0

        for branch, polarity in self.input_devices:
            isum += polarity * branch.state

        denom = self.I * self.R - self.state - self.R * isum

        if denom != 0.0:

            delta_time = self.C * self.R * self.granularity / denom

            if delta_time > 0.0:
                self.tnext = time + delta_time
                self.sign = 1.0

            elif delta_time < 0.0:
                self.tnext = time + -delta_time
                self.sign = -1.0

            else:
                self.tnext = time + 1e-15


class LimGround(qdevs.QdevsDevice):

    def __init__(self):

        qdevs.QdevsDevice.__init__(self)



class LimBranch(qdevs.QdevsDevice):

    def __init__(self, L, R, V, granularity=1e-3, i0=0.0):

        qdevs.QdevsDevice.__init__(self, granularity, i0)

        self.R = R
        self.L = L
        self.V = V

        self.nodei = None
        self.nodej = None

    def connect_nodes(self, nodei, nodej):

        self.connect_outputs(nodei, nodej)
        nodei.connect_outputs(self)
        nodej.connect_outputs(self)

        self.nodei = nodei
        self.nodej = nodej

        nodei.connect_branch(self, polarity=1.0)
        nodej.connect_branch(self, polarity=-1.0)

    def evaluate(self, time):

        self.tnext = float("inf")

        vij = self.nodei.state - self.nodej.state

        denom = self.V + vij - self.R * self.i

        if denom != 0.0:

            delta_time = self.L * self.granularity / denom

            if delta_time > 0.0:
                self.tnext = time + delta_time
                self.trajectory_direction = 1.0

            elif delta_time < 0.0:
                self.tnext = time + -delta_time
                self.trajectory_direction = -1.0

            else:
                self.tnext = time + 1e-15

