
import numpy as np
from scipy.signal import resample
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

import qdevs
import qdevices


def resample(times, values, tf, npoints=500):

    values.append(values[-1])
    times.append(tf)
    f = interp1d(times, values, kind='zero')
    times2 = np.linspace(times[0], times[-1], npoints)
    values2 = f(times2)
    return times2, values2


def test1():

    sys = qdevs.QdevsSystem()

    sq = qdevices.SquareWaveSource(x1=1.0, x2=0.0, t1=2.0, t2=2.0)
    intg = qdevices.Integrator(gain=1.0, granularity=0.1, x0=0.0)
    ode = qdevices.DifferetialEquation(a2=1.0, a1=1.0, a0=1.0, granularity=0.1, x0=0.0)

    sys.add_devices(sq, intg, ode)

    sq.connect_outputs(intg, ode)

    tf = 10.0

    sys.initialize()
    sys.run(tf)

    plt.figure()

    x, t = qdevs.resample(sq.time_history, sq.state_history, tf)
    plt.plot(x, t, label="Sq Wave Output")

    x, t = qdevs.resample(intg.time_history, intg.state_history, tf)
    plt.plot(x, t, label="Integrator Output")

    x, t = qdevs.resample(ode.time_history, ode.state_history, tf)
    plt.plot(x, t, label="ODE Output")

    plt.legend()
    plt.grid()
    plt.show()
    


def test2():

    sys = QdevsSystem()


    dv = 0.01
    di = 0.01

    ground = LimGround()

    node1 = LimNode(C=0.1, R=1.0e-1, I=1.0, v0=0.0, granularity=dv)
    node2 = LimNode(C=1.0, R=1.0e-1, I=1.0, v0=0.0, granularity=dv)
    node3 = LimNode(C=10.0, R=1.0e-1, I=1.0, v0=0.0, granularity=dv)

    branch1 = LimBranch(L=0.1, R=1.0e-1, V=1.0, i0=0.0, granularity=di)
    branch2 = LimBranch(L=1.0, R=1.0e-1, V=1.0, i0=0.0, granularity=di)
    branch3 = LimBranch(L=10.0, R=1.0e-1, V=1.0, i0=0.0, granularity=di)
    branch4 = LimBranch(L=10.0, R=1.0e-1, V=1.0, i0=0.0, granularity=di)

    sys.add_devices(node1, node2, node3, branch1, branch2, branch3,
                    branch4, ground)

    branch1.connect(ground, node1)
    branch2.connect(node1, node2)
    branch3.connect(node1, node3)
    branch4.connect(node2, node3)

    sys.initialize()
    sys.run(500)
    node2.R = 0.01
    sys.run(1000)

    plt.figure()
    plt.plot(*resample_output(node1.times, node1.voltages, 1000, 1000), label="v1")
    plt.plot(*resample_output(node2.times, node2.voltages, 1000, 1000), label="v2")
    plt.plot(*resample_output(node3.times, node3.voltages, 1000, 1000), label="v2")
    plt.legend()
    plt.show()


def test3():

    sys = QdevsSystem()

    grid_size = 20

    nodes = []
    for i in range(grid_size*grid_size):
        if i == 0:
            nodes.append(LimGround())
            continue
        node = LimNode(C=1.0, R=1.0, I=0.0, v0=0.0, granularity=0.1)
        nodes.append(node)

    branches = {}
    for j in range(grid_size):
        for i in range(grid_size):
            if i == 0 and j == 0:
                branch1 = LimBranch(L=1.0, R=1.0, V=100.0, i0=0.0, granularity=0.1)
            else:
                branch1 = LimBranch(L=1.0, R=1.0, V=10.0, i0=0.0, granularity=0.1)
            branch2 = LimBranch(L=1.0, R=1.0, V=50.0, i0=0.0, granularity=0.1)
            branches[(i*j, i*j+1)] = branch1
            branches[(i*j, i*j+grid_size)] = branch2

    sys.add_devices(*nodes)
    sys.add_devices(*branches.values())

    for (i, j), branch in branches.items():
        branch.connect(nodes[i], nodes[j])

    sys.initialize()
    sys.run(10)

    plt.figure()
    for i in range(grid_size)[1:]:
        plt.plot(*resample_output(nodes[grid_size*i].times,
                  nodes[grid_size*i].voltages, 500, 10),
                  label="v({0},{0})".format(i))
    plt.legend()
    plt.show()


if __name__ == "__main__":

    test1()
    #test2()
    #test3()
