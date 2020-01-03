
from math import pi, sin, cos
import liqss
from matplotlib import pyplot as plt


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

    # machine parameters:

    f = 50.0
    omega = 2*pi*f
    Ld = 7.0e-3
    Ll = 2.5067e-3
    Lm = 6.6659e-3
    LD = 8.7419e-3
    LF = 7.3835e-3
    Lq = 5.61e-3
    MQ = 4.7704e-3
    Rs = 1.6e-3
    RF = 9.845e-4
    RD = 0.11558
    RQ = 0.0204
    n = 1.0
    J = 2.812e4
    Cr0 = -2.65e5
    Vfd = 20.0e3
    VF = 10.3

    # simulation parameters:

    dqmin = 1e-7
    dqmax = 1e-2
    dqerr = 0.005
    sr = 80
    tstop = 10.0

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
    b22 = (Lq + Ll) / det2

    # intial conditions:

    fdr0 = 63.662
    fqr0 = 0.0
    fF0 = 72.9852
    fD0 = 65.4777
    fQ0 = 8.8528e-06
    ang0 = 0.0
    wm0 = 314.159

    # derivative functions:

    def dfdr():
       return Vfd * sin(ang.q) - Rs * (a11 * fdr.q + a12 * fF.q + a13 * fD.q) + wm.q * fqr.q

    def dfqr():
        return Vfd * cos(ang.q) - Rs * (b11 * fqr.q + b12 * fQ.q) - wm.q * fdr.q

    def dfF():
        return VF - (a21 * fdr.q + a22 * fF.q + a23 * fD.q) * RF

    def dfD():
        return -(a31 * fdr.q + a32 * fF.q + a33 * fD.q) * RD

    def dfQ():
        return -(b21 * fqr.q + b22 * fQ.q) * RQ

    def dwm():
        return n/J * ((b11 * fqr.q + b12 * fQ.q) * fdr.q - (a11 * fdr.q + a12 * fF.q + a13 * fD.q) * fqr.q - cr.q)

    def dang():
        return wm.q - omega

    gen = liqss.Module("genset", print_time=True)

    cr = liqss.Atom("cr", source_type=liqss.SourceType.RAMP, x1=0.0, x2=Cr0, t1=0.1, t2=5.1, dq=1e4, units="N.m")
    fdr = liqss.Atom("fdr", x0=fdr0, func=dfdr, dmax=sr, units="Wb", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fqr = liqss.Atom("fqr", x0=fqr0, func=dfqr, dmax=sr, units="Wb", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fF = liqss.Atom("fF", x0=fF0, func=dfF, units="Wb", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fD = liqss.Atom("fD", x0=fD0, func=dfD, units="Wb", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    fQ = liqss.Atom("fQ", x0=fQ0, func=dfQ, units="Wb", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    wm = liqss.Atom("wm", x0=wm0, func=dwm, units="rad/s", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)
    ang = liqss.Atom("ang", x0=ang0, func=dang, units="rad", dqmin=dqmin, dqmax=dqmax, dqerr=dqerr)

    fdr.connects(fdr, fqr, fF, fD, wm, ang)
    fqr.connects(fdr, fqr, fF, fD, wm, ang)
    fF.connects(fdr, fqr, fF, fD, wm, ang)
    fD.connects(fdr, fqr, fF, fD, wm, ang)
    fQ.connects(fdr, fqr, fF, fD, wm, ang)

    wm.connects(fqr, fdr, fF, fD, fQ, cr)
    ang.connects(wm)

    gen.add_atoms(cr, fdr, fqr, fF, fD, fQ, wm, ang)

    # simulation:

    gen.initialize()
    gen.run_to(tstop)
    #gen.run_to(tstop, fixed_dt=1e-4)

    r = 4
    c = 2
    i = 0

    plt.figure()

    for i, (name, atom) in enumerate(gen.atoms.items()):

        plt.subplot(r, c, i+1)
        plt.plot(atom.tzoh, atom.qzoh, 'b-')
        #plt.plot(atom.tout, atom.qout, 'k.')
        plt.xlabel("t (s)")
        plt.ylabel(name + " (" + atom.units + ")")

    plt.show()


if __name__ == "__main__":

    #simple()
    #delegates()
    genset()