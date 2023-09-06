
import matplotlib.pyplot as plt
import numpy as np

filename = "sheardata.txt"

data = np.loadtxt(filename)

x = data[:, 0]
y = data[:, 1]
shear_stress = data[:, 2]


indices_to_plot = np.where((x + y) % 10 == 0)

x_filtered = x[indices_to_plot]
y_filtered = y[indices_to_plot]
shear_stress_filtered = shear_stress[indices_to_plot]

plt.figure(figsize=(8, 6))

scatter = plt.scatter(x_filtered, y_filtered, c=shear_stress_filtered, cmap='viridis', s=10)
plt.colorbar(label='Shear Stress')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Shear Stress Distribution')



plt.show()
