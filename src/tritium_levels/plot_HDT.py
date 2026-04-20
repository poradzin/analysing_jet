import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("outputs/tritium_levels.txt")

pulse    = data[:, 0].astype(int)
tritium  = data[:, 1]
deuterium = data[:, 2]
hydrogen = data[:, 3]

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(pulse, tritium, marker='o')
ax1.set_ylabel('Max tritium')
ax1.grid(True)

ax2.plot(pulse, hydrogen, marker='o', color='orange')
ax2.set_ylabel('Hydrogen')
ax2.set_xlabel('Pulse')
ax2.grid(True)

plt.tight_layout()
plt.show()
