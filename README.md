# Python QDEVS Library
Quantized Discrete Event Simulation (QDEVS) engine with Latency Inserion Method (LIM) electrical circuit tolopogy support.

## Example

![plot1.png](img/plot1.png?raw=true "RLC Circuit Simulation Output")

```python

from matplotlib import pyplot as plt
from qdevs import *
from qdevslim import *

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

```
