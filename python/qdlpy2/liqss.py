
from math import sin, asin, pi
from matplotlib import pyplot as plt
from collections import OrderedDict as odict


_EPS = 1.0e-9
_INF = float('inf')
_MAXITER = 1000


class SourceType:

    NONE = "NONE"
    CONSTANT = "CONSTANT"
    STEP = "STEP"
    SINE = "SINE"
    PWM = "PWM"
    RAMP = "RAMP"
    FUNCTION = "FUNCTION"


class Atom(object):

    def __init__(self, name, a=0.0, b=0.0, c=0.0, source_type=SourceType.NONE,
                 x0=0.0, x1=0.0, x2=0.0, xa=0.0, freq=0.0, phi=0.0, func=None,
                 duty=0.0, t1=0.0, t2=0.0, dq=None, dqmin=None, dqmax=None,
                 dqerr=None, dtmin=None, dmax=1e5, units="", srcfunc=None,
                 srcdt=None, output_scale=1.0):

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

        # cache snie wave ta function params:

        self.omega = 2.0 * pi * self.freq

        if self.freq:
            self.T = 1.0 / self.freq

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

        # playback function:
        self.srcfunc = srcfunc
        self.srcdt = srcdt

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

        self.implicit = True

        self.output_scale = output_scale

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

        if self.source_type == SourceType.FUNCTION:
            self.x = self.srcfunc()
            self.q = self.x
            self.q0 = self.x
        else:
            self.x = self.x0
            self.q = self.x0
            self.q0 = self.x0

        if self.dqmin is None:
            self.dqmin = self.sys.dqmin

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
        
        if self.source_type == SourceType.CONSTANT:

            self.x = self.x0

        elif self.source_type == SourceType.FUNCTION:

            self.x = self.srcfunc()

        elif self.source_type == SourceType.NONE:

            self.x += self.d * (self.time - self.tlast)

        elif self.source_type == SourceType.RAMP:

            if self.time <= self.t1:
                self.x = self.x1
            elif self.time <= self.t2:
                self.x = self.x1 + (self.time - self.t1) * self.d 
            else:
                self.x = self.x2

        elif self.source_type == SourceType.STEP:

            if self.time < self.t1:
                self.x = self.x0
            else:
                self.x = self.x1

        elif self.source_type == SourceType.SINE:

            if self.time >= self.t1:
                self.x = self.x0 + self.xa * sin(self.omega * self.time + self.phi)
            else:
                self.x = self.x0

        self.tlast = self.time

    def quantize(self):
        
        interp = False
        change = False

        self.d0 = self.d

        if self.source_type in (SourceType.FUNCTION, SourceType.STEP, SourceType.SINE):

            self.q = self.x

        elif self.source_type in (SourceType.NONE, SourceType.RAMP):

            if self.x >= self.qhi:

                self.q = self.qhi
                self.qlo += self.dq
                change = True

            elif self.x <= self.qlo:

                self.q = self.qlo
                self.qlo -= self.dq
                change = True

            self.qhi = self.qlo + 2.0 * self.dq

            if change and self.implicit:  # we've ventured out of (qlo, qhi) bounds

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

        if self.source_type == SourceType.FUNCTION:

            self.tnext = self.time + self.srcdt

        elif self.source_type == SourceType.NONE:

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

        elif self.source_type == SourceType.STEP:

            if self.time < self.t1:
                self.tnext = self.t1
            else:
                self.tnext = _INF

        elif self.source_type == SourceType.SINE:

            if self.time < self.t1:

                self.tnext = self.t1

            else: 

                w = self.time % self.T             # cycle time
                t0 = self.time - w                 # cycle start time
                theta = self.omega * w + self.phi  # wrapped angular position

                # value at current time w/o dc offset:
                x = self.xa * sin(2.0 * pi * self.freq * self.time)

                # determine next transition time. Saturate at +/- xa:
            
                if theta < pi/2.0:  # quadrant I
                    self.tnext = t0 + (asin(min(1.0, (x + self.dq)/self.xa))) / self.omega

                elif theta < pi:  # quadrant II
                    self.tnext = t0 + self.T/2.0 - (asin(max(0.0, (x - self.dq)/self.xa))) / self.omega

                elif theta < 3.0*pi/2:  # quadrant III
                    self.tnext = t0 + self.T/2.0 - (asin(max(-1.0, (x - self.dq)/self.xa))) / self.omega

                else:  # quadrant IV
                    self.tnext = t0 + self.T + (asin(min(0.0, (x + self.dq)/self.xa))) / self.omega

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


class ComplexAtom(Atom):

    def __init__(self, *args, freq1=1.0, **kwargs):

        Atom.__init__(self, *args, **kwargs)

        self.x0 = complex(self.x0, 0.0)

        # other member variables:
        self.qlo = complex(0.0, 0.0)   
        self.qhi = complex(0.0, 0.0)     
        self.x = self.x0    
        self.d = complex(0.0, 0.0)       
        self.d0 = complex(0.0, 0.0)      
        self.q = self.x0      
        self.q0 = self.x0 
        
        self.freq1 = freq1
        
    def initialize(self, t0):

        self.tlast = t0
        self.time = t0
        self.tnext = _INF

        self.x = self.x0
        self.q = self.x0
        self.q0 = self.x0
        self.dq = self.dqmin
        self.qhi = complex(self.q.real + self.dq, self.q.imag + self.dq)
        self.qlo = complex(self.q.real - self.dq, self.q.imag - self.dq)
        
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

        self.dint()
        self.quantize()

        self.ta()

        # trigger external update if quantized output changed:
        
        if self.q != self.q0:
            self.save()
            self.q0 = self.q
            self.trigger()
            self.update_dq()

    def dint(self):
        
        self.x = complex(self.x.real + self.d.real * (self.time - self.tlast), 
                         self.x.imag + self.d.imag * (self.time - self.tlast))

        self.tlast = self.time

    def step(self, time):

        self.time = time
        self.updates += 1
        self.d = self.f(self.x)
        self.dint()
        self.q = self.x
        self.save()
        self.q0 = self.q

    def quantize(self):
        
        interp = False
        change = False

        self.d0 = self.d

        if self.x.real >= self.qhi.real:

            self.q = complex(self.qhi.real, self.q.imag)
            self.qlo = complex(self.qlo.real + self.dq, self.qlo.imag) 
            change = True

        elif self.x.real <= self.qlo.real:

            self.q = complex(self.qlo.real, self.q.imag)
            self.qlo = complex(self.qlo.real - self.dq, self.qlo.imag) 
            change = True

        if self.x.imag >= self.qhi.imag:

            self.q = complex(self.q.real, self.qhi.imag)
            self.qlo = complex(self.qlo.real, self.qlo.imag + self.dq) 
            change = True

        elif self.x.imag <= self.qlo.imag:

            self.q = complex(self.q.real, self.qlo.imag)
            self.qlo = complex(self.qlo.real, self.qlo.imag - self.dq)
            change = True

        self.qhi = complex(self.qlo.real + 2.0 * self.dq, self.qlo.imag + 2.0 * self.dq)

        if change and self.implicit:  # we've ventured out of (qlo, qhi) bounds

            self.d = self.f(self.q)

            # if the derivative has changed signs, then we know 
            # we are in a potential oscillating situation, so
            # we will set the q such that the derivative ~= 0:

            if (self.d.real * self.d0.real) < 0 or (self.d.imag * self.d0.imag) < 0: 
                flo = self.f(self.qlo) 
                fhi = self.f(self.qhi)
                if flo != fhi:
                    a = (2.0 * self.dq) / (fhi - flo)
                    self.q = self.qhi - a * fhi
                    interp = True

        return interp

    def ta(self):

        if self.d.real > 0.0:
            treal = self.time + (self.qhi.real - self.x.real) / self.d.real
        elif self.d.real < 0.0:
            treal = self.time + (self.qlo.real - self.x.real) / self.d.real
        else:
            treal = _INF

        if self.d.imag > 0.0:
            timag = self.time + (self.qhi.imag - self.x.imag) / self.d.imag
        elif self.d.imag < 0.0:
            timag = self.time + (self.qlo.imag - self.x.imag) / self.d.imag
        else:
            timag = _INF

        self.tnext = min(treal, timag)
        self.tnext = max(self.tnext, self.tlast + self.dtmin)

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
            
        self.dq = min(self.dqmax, max(self.dqmin, abs(self.dqerr * abs(self.q)))) 
            
        self.qlo = self.q - complex(self.dq, self.dq) 

        self.qhi = self.q + complex(self.dq, self.dq) 


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

    def run_to(self, tstop, fixed_dt=None, verbose=False):

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
            if verbose and self.time-last_print_time > 0.01:
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


