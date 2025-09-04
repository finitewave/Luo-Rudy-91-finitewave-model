import numpy as np
import matplotlib.pyplot as plt

from implementation.luo_rudy_91_0d import LuoRudy910D, Stimulation


stimulations = [Stimulation(t_start=0.1, duration=1, amplitude=100.0)]
t_max = 500.0

model = LuoRudy910D(dt=0.01, stimulations=stimulations)
model.run(t_max=t_max)

time = np.arange(0, t_max, model.dt)
plt.plot(time, model.history['u'])
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (u)')
plt.title('0D Luo-Rudy 91 Simulation')
plt.grid()
plt.show()

