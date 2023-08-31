import matplotlib.pyplot as plt
import math

# Read data from a text file (assuming two columns: x-values and y-values)
data_file = 'comparison.txt'
x_values = []
y_values = []
z_values = []

with open(data_file, 'r') as file:
    for line in file:
        x, ux ,uy,uexact= map(float, line.split())
        x_values.append(x)
        y_values.append(math.sqrt(ux*ux + uy*uy))
        z_values.append(uexact)

# Plot the data
plt.figure(figsize=(8, 6))  # Set the figure size

plt.plot(x_values, z_values, marker='o', linestyle='--', color='yellow', label='Exact Profile')

plt.plot(x_values, y_values, marker='x', color='black', label='Lattice Boltzman Profile', linestyle=None,linewidth = 0 )

plt.xlabel('y/H')  # X-axis label
plt.ylabel('Velocity')  # Y-axis label
plt.title('2D Velocity Plot for Poiseulle Flow')  # Plot title
plt.legend()  # Show legend

# Save the plot as a PNG image
plt.savefig('plot.png')

# Show the plot
plt.show()
