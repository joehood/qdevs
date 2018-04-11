# qdevs
Quantized Discrete Event Simulator witten in Python

## Example code

```python
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

```
