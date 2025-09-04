import numpy as np
import matplotlib.pyplot as plt

from implementation import Model0D, Stimulation


stimulations = [Stimulation(t_start=0.1, duration=0.2, amplitude=1.0)]
t_max = 100.0

model = Model0D(dt=0.01, stimulations=stimulations)
model.run(t_max=t_max)

time = np.arange(0, t_max, model.dt)
plt.plot(time, model.history['u'])
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (u)')
plt.title('0D Model Simulation')
plt.grid()
plt.show()

