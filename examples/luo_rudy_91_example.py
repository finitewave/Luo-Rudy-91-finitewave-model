import numpy as np
import matplotlib.pyplot as plt

from implementation.luo_rudy_91_0d import LuoRudy910D, Stimulation


stimulations = [Stimulation(t_start=100., duration=1., amplitude=50.0)]
t_max = 600.0

model = LuoRudy910D(dt=0.01, stimulations=stimulations)
model.run(t_max=t_max)

# fig = plt.figure(figsize=(10, 5))
plt.plot(model.times, model.history['u'], lw=2)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (u)')
plt.title('0D Luo-Rudy 91 Simulation')
plt.grid()
plt.show()

# fig.savefig('luo_rudy_91_ap.png', dpi=300)
