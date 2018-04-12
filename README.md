# Python QDEVS Library
Quantized Discrete Event Simulation (QDEVS) engine with Latency Inserion Method (LIM) electrical circuit tolopogy support.

## Simple Example

```python
from matplotlib import pyplot as plt
from qdevs import *

sys = QdevsSystem(0.1)

sq = SquareWaveSource(x1=1.0, x2=-1.0, t1=2.0, t2=2.0)
intg = Integrator(gain=1.0)
ode = DifferentialEquation(a=-1.0, b=1.0)

sys.add_devices(sq, intg, ode)
sq.connect_outputs(intg, ode)

tf = 6.0
sys.initialize()
sys.run(tf)

plt.figure()
plt.plot(sq.time_history, sq.state_history)
plt.plot(intg.time_history, intg.state_history)
plt.plot(ode.time_history, ode.state_history)
plt.show()

```

## LIM Circuit Simulation Example 

```python
sys = QdevsLimSystem(0.001, 0.001);

n1 = sys.add_node(C=1.0, R=1.0, I=1.0)
g1 = sys.add_ground()
b1 = sys.add_branch(g1, n1, L=1.0, R=1.0, V=1.0)

tf = 10.0

sys.initialize()
sys.run(tf)

plt.figure()
plt.plot(n1.time_history, n1.state_history)
plt.plot(b1.time_history, b1.state_history)
plt.show()

```
