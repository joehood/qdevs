"""Test simulations for qdevs."""

from matplotlib import pyplot as plt

from qdevs import *
from qdevslim import *


def test1():

    sys = QdevsSystem(granularity=0.1)

    src = SquareWaveSource(1, -1, 2, 2)
    intg = Integrator(gain=1.0)
    ode = DifferentialEquation(a=-1.0, b=1.0)

    sys.add_devices(src, intg, ode)
    src.connect_outputs(intg, ode)

    tf = 6.0
    sys.initialize()
    sys.run(tf)

    plt.figure()

    plt.plot(*resample(src.time_history, src.state_history, tf), label="Source Output")
    plt.plot(src.time_history, src.state_history, "k.")

    plt.plot(*resample(intg.time_history, intg.state_history, tf), label="Integrator Output")
    plt.plot(intg.time_history, intg.state_history, "k.")

    plt.plot(*resample(ode.time_history, ode.state_history, tf), label="ODE Output")
    plt.plot(ode.time_history, ode.state_history, "k.")

    plt.legend()
    plt.grid()
    plt.show()
    

def test2():

    sys = QdevsLimSystem(voltage_granularity=0.05, current_granularity=0.05);

    n1 = sys.add_node(C=1.0, R=1.0, I=1.0)
    g1 = sys.add_ground()
    b1 = sys.add_branch(nodei=g1, nodej=n1, L=1.0, R=1.0, V=1.0)
    
    tf = 10.0

    sys.initialize()
    sys.run(tf)

    dt = 0.03
    sys.ss.initialize(dt)
    t, y = sys.ss.run(tf)

    plt.figure()
    plt.plot(*resample(n1.time_history, n1.state_history, tf), label="v (V)")
    plt.plot(t, y[0], 'c--')
    plt.plot(n1.time_history, n1.state_history, "k.")
    plt.plot(*resample(b1.time_history, b1.state_history, tf), label="i (A)")
    plt.plot(t, y[1], 'm--')
    plt.plot(b1.time_history, b1.state_history, "k.")
    plt.legend()
    plt.show()


def test3():

    sys = QdevsLimSystem(0.1, 0.1);

    g1 = sys.add_ground()

    n1 = sys.add_node(C=0.1, R=1.0e-1, I=1.0)
    n2 = sys.add_node(C=1.0, R=1.0e-1, I=1.0)
    n3 = sys.add_node(C=10.0, R=1.0e-1, I=1.0)

    b1 = sys.add_branch(g1, n1, L=0.1, R=1.0e-1, V=1.0)
    b2 = sys.add_branch(n1, n2, L=1.0, R=1.0e-1, V=1.0)
    b3 = sys.add_branch(n1, n3, L=10.0, R=1.0e-1, V=1.0)
    b4 = sys.add_branch(n2, n3, L=10.0, R=1.0e-1, V=1.0)
    
    tf = 15.0
    sys.initialize()
    sys.run(tf)

    #n2.R = 0.01

    #tf = 1000.0
    #sys.run(tf)

    plt.figure()
    plt.plot(*resample(n1.time_history, n1.state_history, tf), label="v1 (V)")
    plt.plot(n1.time_history, n1.state_history, "k.")
    plt.plot(*resample(n2.time_history, n2.state_history, tf), label="v2 (V)")
    plt.plot(n2.time_history, n2.state_history, "k.")
    plt.plot(*resample(n3.time_history, n3.state_history, tf), label="v3 (V)")
    plt.plot(n3.time_history, n3.state_history, "k.")
    plt.legend()
    plt.show()


def test4():

    a = np.matrix([[-1.0, 1.0], [-1.0, -1.0]])
    b = np.matrix([[1.0, 0.0], [0.0, 1.0]])
    u = np.matrix([[1.0], [1.0]])

    sys = StateSpace(a, b, u0=u)

    dt = 0.1
    tf = 10.0

    sys.initialize(dt)
    t, y = sys.run(tf)

    plt.plot(t, y[0])
    plt.plot(t, y[1])
    plt.show()


if __name__ == "__main__":

    #test1()
    test2()
    #test3()
    #test4()
