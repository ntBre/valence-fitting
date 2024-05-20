import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("natoms.dat", dtype=int)
small_data = [d for d in data if d < 100]


fig, ax1 = plt.subplots()

plt.title("Fragment size distribution")
plt.xlabel("Number of atoms")
plt.ylabel("Count")

left, bottom, width, height = [0.4, 0.4, 0.4, 0.3]
ax2 = fig.add_axes([left, bottom, width, height])

ax1.hist(data)
ax2.hist(small_data)

ax2.title.set_text("<100 atoms")

plt.savefig("figs/natoms.png", dpi=300)
