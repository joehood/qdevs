
from matplotlib import pyplot as plt
from collections import OrderedDict as odict

_EPS = 1.0e-9
_INF = float('inf')
_MAXITER = 1000


class SourceType:

    NONE = "NONE"
    CONSTANT = "CONSTANT"
    SINE = "SINE"
    PWM = "PWM"
    RAMP = "RAMP"


class Atom(object):

    def __init__(self, name, a=0.0, b=0.0, c=0.0, source_type=SourceType.NONE,
                 x0=0.0, x1=0.0, x2=0.0, xa=0.0, freq=0.0, phi=0.0, func=None, 
                 duty=0.0, t1=0.0, t2=0.0, dq=None, dqmin=None, dqmax=None,
                 dqerr=None, dtmin=None, dmax=1e5, units=""):

        self.name = name

        self.a = a
        self.b = b
        self.c = c

        self.source_type = source_type
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2
        self.xa = xa
        self.freq = freq
        self.phi = phi
        self.duty = duty
        self.t1 = t1
        self.t2 = t2

        self.ainv = 0.0
        if self.a > 0:
            self.ainv = 1.0 / a

        self.ramp_slope = 0.0
        if (self.t2 - self.t1) > 0:
            self.ramp_slope = (self.x2 - self.x1) / (self.t2 - self.t1)

        self.recieve_from = {}
        self.broadcast_to = []

        self.sys = None  # parent LiqssSystem reference

        # limits: 
        if dq:
            self.dq = dq  
            self.dqmin = dq
            self.dqmax = dq
            self.dqerr = 0.0
        else:
            self.dq = None

        self.dmax = dmax

        if dqmin:
            self.dqmin = dqmin
        elif not dq:
            self.dqmin = None

        if dqmax:
            self.dqmax = dqmax
        elif not dq:
            self.dqmax = None

        if dqerr:
            self.dqerr = dqerr
        elif not dq:
            self.dqerr = None

        if dtmin:
            self.dtmin = dtmin
        else:
            self.dtmin = None

        # delegate derivative function:
        self.func = None
        if func:
            self.func = func

        # other member variables:
        self.qlo = 0.0   
        self.qhi = 0.0    
        self.tlast = 0.0  
        self.tnext = 0.0  
        self.x = x0    
        self.d = 0.0      
        self.d0 = 0.0     
        self.q = x0      
        self.q0 = x0     
        self.triggered = False

        self.tout = None  # output times quantized output
        self.qout = None  # quantized output 

        self.tzoh = None  # zero-order hold output times quantized output
        self.qzoh = None  # zero-order hold quantized output 

        self.updates = 0

        if self.source_type == SourceType.RAMP:
            self.x0 = self.x1

        self.units = units

    def connect(self, atom, coeff=0.0):

        self.recieve_from[atom] = coeff
        atom.broadcast_to.append(self)

    def connects(self, *atoms):

        for atom in atoms:
            self.recieve_from[atom] = 0.0
            atom.broadcast_to.append(self)

    def update_coupling(self, atom, coef):

        if atom in self.recieve_from:
            self.recieve_from[atom] = coeff

    def initialize(self, t0):

        self.tlast = t0
        self.time = t0
        self.tnext = _INF

        self.x = self.x0
        self.q = self.x0
        self.q0 = self.x0
        self.dq = self.dqmin
        self.qhi = self.q + self.dq
        self.qlo = self.q - self.dq
        
        self.tout = [self.time]
        self.qout = [self.q0]
        self.nupd = [0]

        self.tzoh = [self.time]
        self.qzoh = [self.q0]

        self.updates = 0

    def update(self, time):

        self.time = time
        self.updates += 1
        self.triggered = False

        self.d = self.f(self.q)
        self.d = max(self.d, -self.dmax)
        self.d = min(self.d, self.dmax)

        self.dint()
        self.quantize()

        #self.d = self.f(self.q)
        #self.d = max(self.d, -self.dmax)
        #self.d = min(self.d, self.dmax)

        self.ta()

        # trigger external update if quantized output changed:
        
        if self.q != self.q0:
            self.save()
            self.q0 = self.q
            self.trigger()
            self.update_dq()

    def step(self, time):

        self.time = time
        self.updates += 1
        self.d = self.f(self.x)
        self.dint()
        self.q = self.x
        self.save()
        self.q0 = self.q

    def dint(self):
        
        if self.source_type == SourceType.NONE:

            self.x += self.d * (self.time - self.tlast)

        elif self.source_type == SourceType.RAMP:

            if self.time <= self.t1:
                self.x = self.x1
            elif self.time <= self.t2:
                self.x = self.x1 + (self.time - self.t1) * self.d 
            else:
                self.x = self.x2

        self.tlast = self.time

    def quantize(self):
        
        interp = False
        change = False

        self.d0 = self.d

        if self.source_type in (SourceType.NONE, SourceType.RAMP):

            if self.x >= self.qhi:

                self.q = self.qhi
                self.qlo += self.dq
                change = True

            elif self.x <= self.qlo:

                self.q = self.qlo
                self.qlo -= self.dq
                change = True

            self.qhi = self.qlo + 2.0 * self.dq

            if change:  # we've ventured out of (qlo, qhi) bounds

                self.d = self.f(self.q)

                # if the derivative has changed signs, then we know 
                # we are in a potential oscillating situation, so
                # we will set the q such that the derivative ~= 0:

                if (self.d * self.d0) < 0:  # if derivative has changed sign
                    flo = self.f(self.qlo) 
                    fhi = self.f(self.qhi)
                    if flo != fhi:
                        a = (2.0 * self.dq) / (fhi - flo)
                        self.q = self.qhi - a * fhi
                        interp = True

        return interp

    def ta(self):

        if self.source_type == SourceType.NONE:

            if self.d > 0.0:
                self.tnext = self.time + (self.qhi - self.x) / self.d
            elif self.d < 0.0:
                self.tnext = self.time + (self.qlo - self.x) / self.d
            else:
                self.tnext = _INF
        
        elif self.source_type == SourceType.RAMP:

            if self.time < self.t1:
                self.tnext = self.t1

            elif self.time < self.t2:
                if self.d > 0.0:
                    self.tnext = self.time + (self.q + self.dq - self.x) / self.d
                elif self.d < 0.0:
                    self.tnext = self.time + (self.q - self.dq - self.x) / self.d
                else:
                    self.tnext = _INF

            else:
                self.tnext = _INF

        else:
            self.tnext = _INF

        self.tnext = max(self.tnext, self.tlast + self.dtmin)

    def f(self, qval):

        d = 0.0

        if self.func:
            return self.func()

        if self.source_type == SourceType.NONE:

            xsum = 0.0
            for atom, coef in self.recieve_from.items():
                xsum += atom.q * coef
            d = self.ainv * (self.c - qval * self.b - xsum)

        elif self.source_type == SourceType.RAMP:

            d = self.ramp_slope

        else:

            d = 0.0

        return d

    def trigger(self):

        for atom in self.broadcast_to:
            if atom is not self:
                atom.triggered = True

    def update_dq(self):

        if not self.dqerr:
            return
        else:
            if self.dqerr <= 0.0:
                return

        if not (self.dqmin or self.dqmax):
            return

        if (self.dqmax - self.dqmin) < _EPS:
            return
            
        self.dq = min(self.dqmax, max(self.dqmin, abs(self.dqerr * self.q))) 
            
        self.qlo = self.q - self.dq
        self.qhi = self.q + self.dq

    def save(self):
    
        if self.time != self.tout[-1]:

            self.tout.append(self.time)           
            self.qout.append(self.q)
            self.nupd.append(self.updates)

            self.tzoh.append(self.time)           
            self.qzoh.append(self.q0)
            self.tzoh.append(self.time)           
            self.qzoh.append(self.q)


class Module(object):

    def __init__(self, name, dqmin=None, dqmax=None, dqerr=None, dtmin=None,
                 dq=None, print_time=False):
        
        self.name = name

        if dq:
            self.dq = dq
        else:
            self.dq = 1.0e-6

        if dqmin:
            self.dqmin = dqmin
        else:
            self.dqmin = dq

        if dqmax:
            self.dqmax = dqmax
        else:
            self.dqmax = dqmin

        if dqerr:
            self.dqerr = dqerr
        else:
            self.dqerr = 0.0

        if dtmin:
            self.dtmin = dtmin
        else:
            self.dtmin = _EPS

        self.print_time = print_time

        self.atoms = odict()
        
        # simulation variables:
        self.tstop = 0.0  # end simulation time
        self.time = 0.0  # current simulation time
        self.iprint = 0

    def add_atoms(self, *atoms):

        for atom in atoms:

            if not atom.dq:
                atom.dq = self.dq

            if not atom.dqmin:
                atom.dqmin = self.dqmin

            if not atom.dqmax:
                atom.dqmax = self.dqmax

            if not atom.dtmin:
                atom.dtmin = self.dtmin

            atom.sys = self
            self.atoms[atom.name] = atom
            setattr(self, atom.name, atom)

    def initialize(self, t0=0.0):

        self.time = t0

        for atom in self.atoms.values():
            atom.initialize(t0)

    def run_to(self, tstop, fixed_dt=None):

        self.tstop = tstop

        print("Simulation started...")

        if fixed_dt:
            dt = fixed_dt
            while(self.time < self.tstop):
                for atom in self.atoms.values():
                    atom.step(self.time)
                    atom.save()
                self.time += dt
            return

        
        #self.print_percentage(header=True)

        for i in range(1):
            for atom in self.atoms.values():
                atom.update(self.time)
                atom.save()

        i = 0
        while i < _MAXITER:
            triggered = False
            for atom in self.atoms.values():
                if atom.triggered:
                    triggered = True
                    atom.update(self.time)
            if not triggered:
                break
            i += 1

        # main simulation loop:

        tlast = self.time
        last_print_time =  self.time

        while self.time < self.tstop:
            self.advance()
            if self.print_time and self.time-last_print_time > 0.01:
                print("t = {0:5.2f} s".format(self.time))
                last_print_time = self.time
            tlast = self.time

        # force time to tstop and do one update at time = tstop:

        self.time = self.tstop

        for atom in self.atoms.values():
            atom.update(self.time)
            atom.save()

        #self.print_percentage()
        #self.print_percentage(footer=True)

    def advance(self):

        tnext = _INF

        for atom in self.atoms.values():
            tnext = min(atom.tnext, tnext)

        self.time = max(tnext, self.time + _EPS)
        self.time = min(self.time, self.tstop)

        for atom in self.atoms.values():
            if atom.tnext <= self.time or self.time >= self.tstop:
                atom.update(self.time)

        i = 0
        while i < _MAXITER:
            triggered = False
            for atom in self.atoms.values():
                if atom.triggered:
                    triggered = True
                    atom.update(self.time)
            if not triggered:
                break
            i += 1

    def print_percentage(self, header=False, footer=False):
        
        if header:
            print("\n\nPercentage Complete:") #\n+" + "-"*98 + "+")
            return

        if footer:
            print("\nDone.\n\n")
            return

        i = int(self.time / self.tstop * 100.0)

        while self.iprint < i:
            print(str(self.iprint) + "%")
            self.iprint += 1


