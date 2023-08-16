import matplotlib.pyplot as plt

# Read data from a text file (assuming two columns: x-values and y-values)
data_file = 'comparison.txt'
x_values = []
y_values = []

with open(data_file, 'r') as file:
    for line in file:
        x, y ,a,b= map(float, line.split())
        x_values.append(x)
        y_values.append(y)

# Plot the data
plt.figure(figsize=(8, 6))  # Set the figure size

plt.plot(x_values, y_values, marker='o', linestyle='-', color='blue', label='Profile')
plt.xlabel('y/H')  # X-axis label
plt.ylabel('Velocity')  # Y-axis label
plt.title('Plot of Data')  # Plot title
plt.legend()  # Show legend

# Save the plot as a PNG image
plt.savefig('plot.png')

# Show the plot
plt.show()
