
from math import pi, sin, cos, atan2, sqrt
from cmath import rect
import liqss
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def simple():

    sys = liqss.Module("simple")

    node1 = liqss.Atom("node1", 1.0, 1.0, 1.0, dq=1e-2)
    branch1 = liqss.Atom("branch1", 1.0, 1.0, 1.0, dq=1e-2)

    node1.connect(branch1, -1.0)
    branch1.connect(node1, 1.0)

    sys.add_atoms(node1, branch1)

    sys.initialize()
    sys.run_to(10.0)

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(node1.tzoh, node1.qzoh, 'b-')
    plt.plot(node1.tout, node1.qout, 'k.')
    plt.subplot(2, 1, 2)
    plt.plot(branch1.tzoh, branch1.qzoh, 'b-')
    plt.plot(branch1.tout, branch1.qout, 'k.')
    plt.show()


def stiffline():

    dq = 0.01

    sys = liqss.Module("stiffline", dqmin=dq, dqmax=dq)

    node1 = liqss.Atom("node1", 1.0, 1.0, 1.0)
    node2 = liqss.Atom("node2", 1.0e3, 1.0, 1.0)
    node3 = liqss.Atom("node3", 1.0e6, 1.0, 1.0)
    branch1 = liqss.Atom("branch1", 1.0, 1.0, 1.0)
    branch2 = liqss.Atom("branch2", 1.0, 1.0, 1.0)

    node1.connect(branch1, -1.0)
    node2.connect(branch2, 1.0)
    node3.connect(branch2, -1.0)

    branch1.connect(node1, 1.0)
    branch2.connect(node2, -1.0)
    branch2.connect(node3, 1.0)

    sys.add_atoms(node1, node2, node3, branch1, branch2)

    sys.initialize()
    sys.run_to(1.0e5)

    plt.figure()

    plt.subplot(3, 2, 1)
    plt.plot(node1.tzoh, node1.qzoh, 'b-')
    plt.plot(node1.tout, node1.qout, 'k.')

    plt.subplot(3, 2, 2)
    plt.plot(branch1.tzoh, branch1.qzoh, 'b-')
    plt.plot(branch1.tout, branch1.qout, 'k.')

    plt.subplot(3, 2, 3)
    plt.plot(node2.tzoh, node2.qzoh, 'b-')
    plt.plot(node2.tout, node2.qout, 'k.')
    
    plt.subplot(3, 2, 4)
    plt.plot(branch2.tzoh, branch2.qzoh, 'b-')
    plt.plot(branch2.tout, branch2.qout, 'k.')
    
    plt.subplot(3, 2, 5)
    plt.plot(node3.tzoh, node3.qzoh, 'b-')
    plt.plot(node3.tout, node3.qout, 'k.')
    plt.show()


def delegates():

    dqmin = 0.01
    dqmax = 0.01
    dqerr = 0.01

    sys = liqss.Module("delegates", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    def f1():
        return 1.0 - sys.node1.q + sys.branch1.q

    def f2():
        return 1.0 - sys.branch1.q - sys.node1.q

    node1 = liqss.Atom("node1", func=f1)
    branch1 = liqss.Atom("branch1", func=f2)

    node1.connect(branch1)
    branch1.connect(node1)

    sys.add_atoms(node1, branch1)

    sys.initialize()
    sys.run_to(10.0)

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(node1.tzoh, node1.qzoh, 'b-')
    plt.plot(node1.tout, node1.qout, 'k.')
    plt.subplot(2, 1, 2)
    plt.plot(branch1.tzoh, branch1.qzoh, 'b-')
    plt.plot(branch1.tout, branch1.qout, 'k.')
    plt.show()


def genset():

    # parametmrs:

    # machine:
    f = 50.0
    wb = 2*pi*f
    Ld = 7.0e-3
    Ll = 2.5067e-3
    Lm = 6.6659e-3
    LD = 8.7419e-3
    LF = 7.3835e-3
    Lq = 5.61e-3
    MQ = 4.7704e-3
    Ra = 0.001
    Rs = 1.6e-3
    RF = 9.845e-4
    RD = 0.11558
    RQ = 0.0204
    n = 1.0
    J = 2.812e4

    Xd = wb * Ld
    Ra = Rs

    # converter:
    Rf = 0.001
    Lf = 0.01
    Cf = 0.01
    Clim = 0.001
    Rlim = 100.0

    # propulsion system:
    Rp = 100.0
    Lp = 1000.0

    # avr control:
    Ka = 10.0/120.0e3
    Ta = 10.0

    Tm_max = -2.65e5
    vb = 20.0e3
    efd0 = 10.3
    efd_cmd = 9.406

    # simulation parametmrs:

    dqmin = 1e-7
    dqmax = 1e-2
    dqerr = 0.005
    sr = 80
    tstop = 5.0

    # derived:

    det1 = Lm*(Lm**2-LF*Lm)+Lm*(Lm**2-LD*Lm)+(Ll+Ld)*(LD*LF-Lm**2)
    det2 = (Lq+Ll)*MQ-MQ**2
    a11 = (LD*LF-Lm**2) / det1
    a12 = (Lm**2-LD*Lm) / det1
    a13 = (Lm**2-LF*Lm) / det1
    a21 = (Lm**2-LD*Lm) / det1
    a22 = (LD*(Ll+Ld)-Lm**2) / det1
    a23 = (Lm**2-(Ll+Ld)*Lm) / det1
    a31 = (Lm**2-LF*Lm) / det1
    a32 = (Lm**2-(Ll+Ld)*Lm) / det1
    a33 = (LF*(Ll+Ld)-Lm**2) / det1
    b11 = MQ / det2
    b12 = -MQ / det2
    b21 = -MQ / det2
    b22 = (Lq+Ll) / det2

    # initial states:

    tm0 = -265000.0
    fdr0 = 55.16705817631649
    fqr0 = -31.81417681249221
    fF0 = 67.02932640721295
    fD0 = 59.52167785464967
    fQ0 = -18.69803603266313
    wr0 = 314.1592653589792
    theta0 = 0.5231522988502871
    vdc0 = 12359.114942024915
    vt0 = 9207.54123275002
    idc0 = 123.58915423566667
    vf0 = 12358.991351636208
    ip0 = 123.58915493130449

    S = sqrt(3/2) * 2 * sqrt(3) / pi

    # algebraic functions:

    def Sd():
        return S * cos(theta.q)

    def Sq():
        return -S * sin(theta.q) 

    def efd():
        return efd0

    def id():
        return a11 * fdr.q + a12 * fF.q + a13 * fD.q 

    def iq():
        return b11 * fqr.q + b12 * fQ.q

    def iF():
        return a21 * fdr.q + a22 * fF.q + a23 * fD.q

    def iD():
        return a31 * fdr.q + a32 * fF.q + a33 * fD.q

    def iQ():
        return b21 * fqr.q + b22 * fQ.q

    def vt():
        return sqrt(Xd**2 + Ra**2) * sqrt((iq() - Sq() * idc.q)**2 + (id() - Sd() * idc.q)**2)

    def ed():
        return vb * sin(theta.q)

    def eq():
        return vb * cos(theta.q)

    def vdc():
        return (Sd() * (Ra * (id() - Sd() * idc.q) - Xd * (iq() - Sq() * idc.q))
             + Sq() * (Ra * (iq() - Sq() * idc.q) + Xd * (id() - Sd() * idc.q)))

    def didc():
        return 1/Lf * (Vdc.q - idc.q * Rf - vf.q)

    def dvf():
        return 1/Cf * (idc.q - ip.q)

    def dip():
        return 1/Lp * (vf.q - ip.q * Rp)

    # derivative functions:

    def dfdr():
       return ed() - Rs * id() + wr.q * fqr.q

    def dfqr():
        return eq() - Rs * iq() - wr.q * fdr.q

    def dfF():
        return efd() - iF() * RF

    def dfD():
        return -iD() * RD

    def dfQ():
        return -iQ() * RQ

    def dwr():
        return (n/J) * (iq() * fdr.q - id() * fqr.q - tm.q)

    def dtheta():
        return wr.q - wb

    def davr():
        #return (1/Ta) * (Ka * sqrt(vd.q**2 + vq.q**2) - avr.q)
        return (1/Ta) * (Ka * vdc.q - avr.q)   #  v = i'*L + i*R    i' = (R/L)*(v/R - i)

    ship = liqss.Module("genset", print_time=True)

    # machine:
    tm    = liqss.Atom("tm", source_type=liqss.SourceType.RAMP, x1=Tm_max, x2=Tm_max*0.8, t1=30.0, t2=60.0, dq=1e4, units="N.m")
    #tm    = liqss.Atom("tm", source_type=liqss.SourceType.CONSTANT, x0=tm0, units="N.m", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    
    fdr   = liqss.Atom("fdr",   x0=fdr0,   func=dfdr,   units="Wb",    dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fqr   = liqss.Atom("fqr",   x0=fqr0,   func=dfqr,   units="Wb",    dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fF    = liqss.Atom("fF",    x0=fF0,    func=dfF,    units="Wb",    dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fD    = liqss.Atom("fD",    x0=fD0,    func=dfD,    units="Wb",    dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fQ    = liqss.Atom("fQ",    x0=fQ0,    func=dfQ,    units="Wb",    dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    wr    = liqss.Atom("wr",    x0=wr0,    func=dwr,    units="rad/s", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    theta = liqss.Atom("theta", x0=theta0, func=dtheta, units="rad",   dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    #idc   = liqss.Atom("idc",   x0=idc0,   func=didc,   units="A",     dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    #vf    = liqss.Atom("vf",    x0=vf0,    func=dvf,    units="V",     dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    #ip    = liqss.Atom("ip",    x0=ip0,    func=dip,    units="A",     dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    Vdc = liqss.Atom("vdc", source_type=liqss.SourceType.FUNCTION, srcfunc=vdc, srcdt=1e-3, x0=vdc0, units="V", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    Vt  = liqss.Atom("vt", source_type=liqss.SourceType.FUNCTION, srcfunc=vt, srcdt=1e-3, x0=vt0, units="V", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    It  = liqss.Atom("it", source_type=liqss.SourceType.FUNCTION, srcfunc=it, srcdt=1e-3, x0=it0, units="A", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    

    fdr.connects(fdr, fqr, fF, fD, wr, theta)
    fqr.connects(fdr, fqr, fF, fD, wr, theta)
    fF.connects(fdr, fqr, fF, fD, wr, theta)
    fD.connects(fdr, fqr, fF, fD, wr, theta)
    fQ.connects(fdr, fqr, fF, fD, wr, theta)
    wr.connects(fqr, fdr, fF, fD, fQ, tm)
    theta.connects(wr)

    ship.add_atoms(tm, fdr, fqr, fF, fD, fQ, wr, theta)
    ship.add_atoms(Vdc, Vt)
    #ship.add_atoms(idc, vf, ip)

    # simulation:

    ship.initialize()
    ship.run_to(30.0, fixed_dt=1e-3)

    Rp = 80.0
    #tm.x0 = tm0*1.2
    ship.run_to(90.0, fixed_dt=1e-3)

    # plot all results:

    r = 4
    c = 4
    i = 0

    plt.figure()

    f = open(r"c:\temp\initcond.txt", "w")

    j = 0

    for i, (name, atom) in enumerate(ship.atoms.items()):

        plt.subplot(r, c, j+1)
        plt.plot(atom.tzoh, atom.qzoh, 'b-')
        #plt.plot(atom.tout, atom.qout, 'k.')
        plt.xlabel("t (s)")
        plt.ylabel(name + " (" + atom.units + ")")
        j += 1

        f.write("\t{}0 = {}\n".format(name, atom.q))

    f.close()
    plt.show()


def shipsys():
        
    # machine parameters:

    f = 60.0
    Ld = 7.0e-3
    Lq = 5.61e-3
    Rs = 1.6e-3
    P = 2.0
    J = 2.812e4
    Tm0 = -2.65e5
    Efd = 20.0e3
    Cf = 0.001
    Clim = 0.001
    Lf = 0.001
    Rf = 0.001
    Lp = 10.0
    Rp = 1.0

    # simulation parameters:

    dqmin = 1e-7
    dqmax = 1e-2
    dqerr = 0.005
    sr = 80
    tstop = 0.1

    # initial states:

    vd0 = 0.0  
    vq0 = 0.0  
    idc0 = 0.0 
    vdc0 = 0.0  
    ip0 = 0.0 
    id0 = 0.0 
    iq0 = 0.0 
    wm0 = 377.0  
    theta0 = 0.0     

    S = sqrt(3/2) * 2 * sqrt(3) / pi
    wb = 2*pi*f

    def Sd():
        return S * cos(theta.q)

    def Sq():
        return S * sin(theta.q)

    # derivative functions:
 
    def dvd():
        return 1/Clim * (id.q - idc.q * Sd()) 

    def dvq():
        return 1/Clim * (iq.q - idc.q * Sq())
 
    def didc():
        return 1/Lf * (vd.q * Sd() + vq.q * Sq() - idc.q * Rf - vdc.q)
 
    def dvdc():
        return 1/Cf * (idc.q - ip.q) 
 
    def dip():
        return 1/Lp * (vdc.q - ip.q * Rp)
 
    def did():
        return 1/Ld * (Efd * sin(theta.q) - id.q * Rs - wm.q * iq.q * Lq - vd.q)
 
    def diq():
        return 1/Lq * (Efd * cos(theta.q) - iq.q * Rs + wm.q * id.q * Ld - vq.q)
 
    def dwm():
        return P/(2*J) * (3/2*P/2 * (id.q * Ld * iq.q - iq.q * Lq * id.q) - Tm0)
 
    def dtheta():
        return wm.q - wb

    ship = liqss.Module("ship", print_time=True)

    #tm = liqss.Atom("tm", source_type=liqss.SourceType.RAMP, x1=0.0, x2=Tm0, t1=0.1, t2=5.1, dq=1e4, units="N.m")
    #tm = liqss.Atom("tm", source_type=liqss.SourceType.CONSTANT, x1=0.0, units="N.m", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    
    vd    = liqss.Atom("vd",    x0=vd0   , func=dvd   , units="V"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)           
    vq    = liqss.Atom("vq",    x0=vq0   , func=dvq   , units="V"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    idc   = liqss.Atom("idc",   x0=id0   , func=didc  , units="A"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    vdc   = liqss.Atom("vdc",   x0=vd0   , func=dvdc  , units="V"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    ip    = liqss.Atom("ip",    x0=ip0   , func=dip   , units="A"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    id    = liqss.Atom("id",    x0=id0   , func=did   , units="A"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    iq    = liqss.Atom("iq",    x0=iq0   , func=diq   , units="A"    , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    wm    = liqss.Atom("wm",    x0=wm0   , func=dwm   , units="rad/s", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    theta = liqss.Atom("theta", x0=theta0, func=dtheta, units="rad"  , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    vd.connects(id, iq)   
    vq.connects(id, iq)    
    idc.connects(vd, vq, id, iq)   
    vdc.connects(id, iq, ip)    
    ip.connects(vdc, ip)    
    id.connects(vd, vq, id, iq)    
    iq.connects(vd, vq, id, iq)    
    wm.connects(id, iq)    
    theta.connects(wm) 

    ship.add_atoms(vd, vq, idc, vdc, ip, id, iq, wm, theta)

    # simulation:

    # no-load startup:

    ship.initialize()
    #ship.run_to(tstop)
    ship.run_to(tstop, fixed_dt=1e-6)

    r = 5
    c = 2
    i = 0

    plt.figure()

    for i, (name, atom) in enumerate(ship.atoms.items()):

        plt.subplot(r, c, i+1)
        plt.plot(atom.tzoh, atom.qzoh, 'b-')
        #plt.plot(atom.tout, atom.qout, 'k.')
        plt.xlabel("t (s)")
        plt.ylabel(name + " (" + atom.units + ")")

    plt.show()


def shipsys2():
        
    # machine parameters:

    Prate = 555.0e6 # V.A
    Vrate = 24.0e3  # V_LLRMS
    freq = 60.0     # Hz
    P = 2 # number of poles

    wb = 2*pi*freq  # base speed

    vb = 4160.0 # base voltage RMS LL

    Lad = 1.66
    Laq = 1.61
    Lo = 0.15
    Ll = 0.15
    Ra = 0.003
    Lfd = 0.165
    Rfd = 0.0006
    L1d = 0.1713
    R1d = 0.0284
    L1q = 0.7252
    R1q = 0.00619
    L2q = 0.125
    R2q = 0.2368

    J = 2.525 # 27548.0

    # derived:
    Lffd = Lad + Lfd
    Lf1d = Lffd - Lfd
    L11d = Lf1d + L1d
    L11q = Laq + L1q
    L22q = Laq + L2q

    # simulation parameters:
    dqmin = 1e-7
    dqmax = 1e-2
    dqerr = 0.005
    sr = 80
    tstop = 0.1

    # initial states:

    efd0 = 92.95
    fd0 = 0.0
    fq0 = 0.0
    fo0 = 0.0
    ffd0 = 0.0
    f1d0 = 0.0
    f1q0 = 0.0
    f2q0 = 0.0
    wr0 = wb
    theta0 = 0.0

    # algebraic equations:

    def id():
        return -(((Lad*Lf1d-L11d*Lad) * ffd.q + (L11d*Lffd-Lf1d**2) * fd.q + (Lad*Lf1d-Lad*Lffd) * f1d.q)
                 / ((L11d*Lffd-Lf1d**2)*Ll+(L11d*Lad-Lad**2)*Lffd-Lad*Lf1d**2+2*Lad**2*Lf1d-L11d*Lad**2))

    def iq():
        return -(((Laq**2-L11q*L22q)*fq.q + (L11q*Laq-Laq**2)*f2q.q + (L22q*Laq-Laq**2) * f1q.q)
                 / ((Laq**2-L11q*L22q)*Ll-Laq**3+(L22q+L11q)*Laq**2-L11q*L22q*Laq))

    def io():
        return -fo.q / Lo

    def ifd():
        return (((L11d*Ll-Lad**2+L11d*Lad) * ffd.q + (Lad*Lf1d-L11d*Lad) * fd.q + (-Lf1d*Ll-Lad*Lf1d+Lad**2) * f1d.q)
                /((L11d*Lffd-Lf1d**2)*Ll+(L11d*Lad-Lad**2)*Lffd-Lad*Lf1d**2+2*Lad**2*Lf1d-L11d*Lad**2))

    def i1d():
        return -(((Lf1d*Ll+Lad*Lf1d-Lad**2) * ffd.q + (Lad*Lffd-Lad*Lf1d) * fd.q + (-Lffd*Ll-Lad*Lffd+Lad**2) * f1d.q)
                 / ((L11d*Lffd-Lf1d**2)*Ll+(L11d*Lad-Lad**2)*Lffd-Lad*Lf1d**2+2*Lad**2*Lf1d-L11d*Lad**2))

    def i1q():
        return -(((Laq**2-L22q*Laq) * fq.q - (Laq*Ll) * f2q.q + (L22q*Ll-Laq**2+L22q*Laq) * f1q.q)
                 / ((Laq**2-L11q*L22q)*Ll-Laq**3+(L22q+L11q)*Laq**2-L11q*L22q*Laq))

    def i2q():
        return -(((Laq**2-L11q*Laq) * fq.q + (L11q*Ll-Laq**2+L11q*Laq) * f2q.q - (Laq*Ll) * f1q.q)
                 /((Laq**2-L11q*L22q)*Ll-Laq**3+(L22q+L11q)*Laq**2-L11q*L22q*Laq)) 

    def Te():
        return fd.q * iq() + fq.q * id()

    def ed():
        return vb * cos(theta.q)

    def eq():
        return -vb * sin(theta.q)

    def eo():
        return 0

    def efd():
        return efd0

    def tm():
        return 0

    # derivative functions:
 
    def dfd():
        # ed = 1/wb * dfd - fq*wr - Ra * id
        return wb * (ed() + fq.q * wr.q + Ra * id())

    def dfq():
        # eq = 1/wb * dfq + fd*wr - Ra * iq
        return wb * (eq() + fd.q * wr.q + Ra * iq())

    def dfo():
        # eo = 1/wb * dfo - Ra * io
        return wb * (eo() + Ra * io())

    def dffd():
        # efd = 1/wb * dffd - Rfd * ifd
        return wb * (efd() + Rfd * ifd())

    def df1d():
        # 0 = 1/wb * df1d - R1d * i1d
        return wb * R1d * i1d()

    def df1q():
        # 0 = 1/wb * df1q - R1q * i1q
        return wb * R1q * i1q()

    def df2q():
        # 0 = 1/wb * df2q - R2q * i2q
        return wb * R2q * i2q()

    def dwr():
        return P/(2*J) * (Te() - tm())
 
    def dtheta():
        return wr.q - wb 

    ship = liqss.Module("ship", print_time=True)

    #tm = liqss.Atom("tm", source_type=liqss.SourceType.RAMP, x1=0.0, x2=Tm0, t1=0.1, t2=5.1, dq=1e4, units="N.m")

    #tm = liqss.Atom("tm", source_type=liqss.SourceType.CONSTANT, x1=0.0, units="N.m", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
              
    fd     = liqss.Atom("fd"   , x0=fd0   , func=dfd   , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fq     = liqss.Atom("fq"   , x0=fq0   , func=dfq   , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fo     = liqss.Atom("fo"   , x0=fo0   , func=dfo   , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    ffd    = liqss.Atom("ffd"  , x0=ffd0  , func=dffd  , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    f1d    = liqss.Atom("f1d"  , x0=f1d0  , func=df1d  , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    f1q    = liqss.Atom("f1q"  , x0=f1q0  , func=df1q  , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    f2q    = liqss.Atom("f2q"  , x0=f2q0  , func=df2q  , units="Wb"   , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    wr     = liqss.Atom("wr"   , x0=wr0   , func=dwr   , units="rad/s", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    theta  = liqss.Atom("theta", x0=theta0, func=dtheta, units="rad"  , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    fd.connects(fd, f1d, ffd, fq, wr)
    fq.connects(fq, f1q, f2q, fd, wr)
    fo.connects(fo)
    ffd.connects(fd, f1d, ffd)
    f1d.connects(fd, f1d, ffd)
    f1q.connects(fq, f1q, f2q)
    f2q.connects(fq, f1q, f2q)
    wr.connects(fd, fq, f1d, ffd, f1q, f2q, wr)
    theta.connects(wr)

    ship.add_atoms(
      fd,    
      fq,    
      fo,    
      ffd,   
      f1d,   
      f1q,   
      f2q,   
      wr,    
      theta)

    # simulation:

    ship.initialize()
    #ship.run_to(tstop)
    ship.run_to(0.0003, fixed_dt=1e-8)

    r = 5
    c = 2
    i = 0

    plt.figure()

    for i, (name, atom) in enumerate(ship.atoms.items()):

        plt.subplot(r, c, i+1)
        plt.plot(atom.tzoh, atom.qzoh, 'b-')
        #plt.plot(atom.tout, atom.qout, 'k.')
        plt.xlabel("t (s)")
        plt.ylabel(name + " (" + atom.units + ")")

    plt.show()


def gencls():

    """ Simplified Sync Machine Model simulation

    
    """

    # parameters:

    H = 3.0
    Kd = 1.0
    fs = 60.0
    ws = 2*pi*fs


    # intial conditions:

    Tm0 = 5.5


    # odes:

    def dtheta():
        return (wr.q - ws) / ws

    # model and state stoms:

    gencls = liqss.Module("gencls", print_time=True)

    theta  = liqss.Atom("theta", x0=theta0, func=dtheta, units="rad"  , dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    theta.connects()

    ship.add_atoms(theta)

    ship.initialize()

    ship.run_to(0.0003, fixed_dt=1e-8)

    r = 1
    c = 2
    i = 0

    plt.figure()

    for i, (name, atom) in enumerate(ship.atoms.items()):

        plt.subplot(r, c, i+1)
        plt.plot(atom.tzoh, atom.qzoh, 'b-')
        plt.plot(atom.tout, atom.qout, 'k.')
        plt.xlabel("t (s)")
        plt.ylabel(name + " (" + atom.units + ")")

    plt.show()



def modulate(times, values, freq, npoints):

    mags = [abs(x) for x in values]
    phs = [atan2(x.imag, x.real) for x in values]

    fmag = interp1d(times, mags, kind='zero')
    fph = interp1d(times, phs, kind='zero')

    times2 = np.linspace(times[0], times[-1], npoints)
    mags2 = fmag(times2)
    phs2 = fph(times2)

    values2 = [mag * sin(2.0*pi*freq*t + ph) for (t, mag, ph) in zip(times2, mags2, phs2)]

    return times2, values2


def plot_cqss(*catoms):

    nplots = len(catoms)

    plt.figure()

    for iplot, catom in enumerate(catoms):

        ax = plt.subplot(nplots, 1, iplot+1)

        l1 = plt.plot(catom.tzoh, [abs(x) for x in catom.qzoh], 'c-', color='lightblue', linewidth=1)
        l2 = plt.plot(catom.tout, [abs(x) for x in catom.qout], 'b.', color="darkblue", markersize=3)
        tmod, qmod = modulate(catom.tout, catom.qout, catom.freq1, 1000)
        l3 = plt.plot(tmod, qmod, 'b-', color='grey', linewidth=0.5, label="x(t) (modulated)")
        plt.xlabel("t(s)")

        plt.ylabel(catom.name + " (" + catom.units + ")", color='blue')

        #ax2 = ax.twinx()
        #l4 = ax2.plot(catom.tout, catom.nupd, 'c-', color='red', linewidth=1, label="cummulative updates")
        #ax2.set_ylabel("updates", color='red')
        #labels = ["|X(t)| qss (zoh)", "|X(t)| qss", "x(t) (modulated)", "cummulative updates"]
        #plt.legend(handles=[l1[0], l2[0], l3[0], l4[0]], labels=labels, loc='lower right')

        labels = ["|X(t)| qss (zoh)", "|X(t)| qss", "x(t) (modulated)"]
        plt.legend(handles=[l1[0], l2[0], l3[0]], labels=labels, loc='lower right')

        ax.spines['left'].set_color('blue')
        ax.tick_params(axis='y', colors='blue', which='both')
        #ax2.spines['right'].set_color('red')
        #ax2.tick_params(axis='y', colors='red', which='both')

    plt.show()


def dynphasor():

    """
             I
            --->
     .---VVV---UUU---o------.------.     
    +|    R     L    |      |      |      +
  E (~)           C ===  G [ ]  H (~) v   V
    -|               |      |      |      -
     '---------------+------'------'
                    _|_
                     -

    VL = ZL * IL + d*IL/dt * L

    E = I'*L + I*(R + jwL) + V
    I = V'*C + V*(G + 1/jwC) + H

    I' = (1/L) * (E - I*(R + jwL)) - V 
    V' = (1/C) * (I - V*(G + 1/jwC)) - H

    """

    dqmin = 0.03
    dqmax = 0.03
    dqerr = 0.01

    f = 2.0
    omega = 2.0*pi*f
    E = rect(1.0, 0.0)
    H = rect(1.0, 0.0)
    R = 10.0
    L = 1.0
    C = 1.0
    G = 10.0
    jwL = complex(0.0, omega*L)
    jwC = complex(0.0, omega*C)

    sys = liqss.Module("dynphasor", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    def dI():
        return (1/L) * (E - branch1.q*(R + jwL)) - node1.q 

    def dV():
        return (1/C) * (branch1.q - node1.q*(G + 1/jwC)) - H

    branch1 = liqss.ComplexAtom("branch1", units="A", func=dI, dq=0.00002, freq1=f)
    node1 = liqss.ComplexAtom("node1", units="V", func=dV, dq=0.00001, freq1=f)
    
    branch1.connect(node1)
    node1.connect(branch1)

    sys.add_atoms(node1, branch1)

    sys.initialize()
    
    sys.run_to(4, verbose=True)

    E = 2.0 * E 

    sys.run_to(8, verbose=True)

    plot_cqss(node1, branch1)


def genphasor():

    # simulation parameters:

    dqmin = 0.03
    dqmax = 0.03
    dqerr = 0.01

    # per unit bases:

    sbase = 100.0
    fbase = 60.0

    # parameters (per unit):

    Ra    = 0.01
    Tdop  = 6.50
    Tdopp = 0.06
    Tqop  = 0.20
    Tqopp = 0.05
    H     = 4.00
    D     = 0.00
    Xd    = 1.80
    Xq    = 1.75
    Xdp   = 0.60
    Xqp   = 0.80
    Xdpp  = 0.30
    Xqpp  = 0.30
    Xl    = 0.15

    Pmech = 750.0 / sbase

    # derived: 
          
    omega0 = 2.0 * pi * fbase
    G = Ra / (Xdpp**2 + Ra**2)
    B = -Xdpp / (Xdpp**2 + Ra**2)

    # initial conditions:

    delta0 = 0.0
    eqpp0  = 0.0
    edpp0  = 0.0
    eqp0   = 0.0
    edp0   = 0.0

    # algebraic functions:

    def Telec():
       return edpp.q*iq() - eqpp.q*id()

    def id():
        return ((Ra * edpp + Xdpp * eqpp) * omega) / ((Xdpp**2 + Ra**2) * omega0)

    def iq():
        return ((Ra * eqpp - Xdpp * edpp) * omega) / ((Xdpp**2 + Ra**2) * omega0)

    def theta():
        return delta.q

    def Isource():
        return complex(cos(delta.q) * id - sin(delta.q) * iq, 
                       cos(delta.q) * iq + sin(delta.q) * id)

    def efd():
        return 1.0

    # derivative functions:

    def ddelta():
        return omega.q * omega0

    def domega():
        return 1/(2*H) * ((Pmech - D*omega.q)/(1 + omega.q) - Telec())

    def deqpp():
        return (edp.q - edpp.q + (Xqp - Xl) * iq()) / Tqopp

    def dedpp():
        return (eqp.q - eqpp.q - (Xdp - Xl) * id()) / Tdopp

    def deqp():
        return (efd() - (Xd-Xdpp)/(Xdp-Xdpp) * eqp.q + (Xd-Xdp)/(Xdp-Xdpp) * eqpp.q) / Tdop

    def dedp():
        return (-(Xq-Xqpp)/(Xqp-Xqpp) * edp.q + (Xq-Xqp)/(Xqp-Xqpp) * edpp.q) / Tqop

    # system object:

    sys = liqss.Module("genphasor", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    # atoms:

    delta = liqss.Atom("delta", x0=delta0, func=ddelta, units="rad",   dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    omega = liqss.Atom("omega", x0=omega0, func=ddelta, units="rad/s", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    eqpp  = liqss.Atom("eqpp",  x0=eqpp0,  func=ddelta, units="Vpu",   dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    edpp  = liqss.Atom("edpp",  x0=edpp0,  func=ddelta, units="Vpu",   dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    eqp   = liqss.Atom("eqp",   x0=eqp0,   func=ddelta, units="Vpu",   dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    edp   = liqss.Atom("edp",   x0=edp0,   func=ddelta, units="Vpu",   dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    sys.add_atoms(delta, omega, eqpp, edpp, eqp, edp)

    sys.initialize()
    
    sys.run_to(0.1, verbose=True, fixed_dt=1.0e-5)

    #sys.run_to(2.0, verbose=True, fixed_dt=1.0e-4)

    # plot all results:

    r = 3
    c = 2
    i = 0

    plt.figure()

    f = open(r"c:\temp\initcond.txt", "w")

    j = 0

    for i, (name, atom) in enumerate(sys.atoms.items()):

        plt.subplot(r, c, j+1)
        plt.plot(atom.tzoh, atom.qzoh, 'b-')
        #plt.plot(atom.tout, atom.qout, 'k.')
        plt.xlabel("t (s)")
        plt.ylabel(name + " (" + atom.units + ")")
        j += 1

        f.write("\t{}0 = {}\n".format(name, atom.q))

    f.close()
    plt.show()


if __name__ == "__main__":

    #simple()
    #stiffline()
    #delegates()
    #genset()
    #shipsys()
    #shipsys2()
    #dynphasor()
    #gencls()
    genphasor()