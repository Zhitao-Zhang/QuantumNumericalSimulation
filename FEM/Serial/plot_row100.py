import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Read only row 100 (1-indexed) => index 99
with open("trace_zaxis_txx.dat") as f:
    for i, line in enumerate(f):
        if i == 99:
            row = np.array([float(x) for x in line.split()])
            break

x = np.arange(row.size)
plt.figure(figsize=(12, 4))
plt.plot(x, row, lw=0.8)
plt.xlabel("column index")
plt.ylabel("amplitude")
plt.title("trace_zaxis_txx.dat  row 100  ({} cols)".format(row.size))
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("trace_zaxis_txx_row100.png", dpi=150)
print("saved trace_zaxis_txx_row100.png  min=%.3e max=%.3e" % (row.min(), row.max()))
