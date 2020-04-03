
import numpy as np
import numpy.linalg as la


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
