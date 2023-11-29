# BTP-1-Sem--2023-2024-LBM
This Repository contains the code files for Lattice Boltzmann Method for our BTP Project.

# Lattice Boltzmann Method for Flow Past a Cylinder

## Introduction
This project simulates the fluid flow past a cylinder using the Lattice Boltzmann Method (LBM), a computational fluid dynamics technique. The previous part of LBM.cpp involves comparison of Navier Stokes with LBM for flow between two parallel plates.

## Features
- 2D fluid flow simulation around a cylindrical obstacle.
- Visualization of velocity fields and vorticity.
- Customizable parameters for different flow scenarios.

## Dependencies
- C++ Compiler (e.g., g++, clang)
- Python 
- NumPy
- Matplotlib

## Getting Started

To use this simulation, first clone the repository to your local machine:

```bash
git clone https://github.com/Kabeer120/BTP-1-Sem--2023-2024-LBM.git
cd BTP-1-Sem--2023-2024-LBM
```
Ensure you have a C++ compiler installed and compile the simulation code:

```bash
g++ -o flowPastCylinder flowPastCylinder.cpp -std=c++11
```
```bash
./flowPastCylinder
```




## Visualizing
- For Flow Between Parellel Plates, run the python file Plot_LBM.py
```bash
python Plot_LBM.py
```
-For Flow Across Cylinder, use the data file generated in Paraview.


## Code Structure

Code Structure
flowPastCylinder.cpp: The simulation operates on a 2D grid, discretizing the fluid domain. It incorporates the D2Q9 model, which means each grid point interacts with its eight neighbors and itself. The program focuses on capturing the dynamics of fluid flow, particularly the interaction with a cylindrical obstacle, and demonstrates features like initialization, collision processes, boundary conditions, and data generation for visualization.

Key Components
Grid and Obstacle Initialization: The simulation grid is defined with dimensions nx and ny. Additionally, a circular obstacle is positioned within the flow field, characterized by its center (circleX, circleY) and radius (circleR). The grid includes 'ghost' nodes to facilitate boundary condition implementation.

Physical and Simulation Parameters: Key parameters such as the speed of sound (cs), Mach number (Ma), Reynolds number (Re), and relaxation time (tau) are initialized. These parameters are crucial in determining the flow characteristics and the stability of the simulation.

Memory Allocation for Simulation Data: Arrays are dynamically allocated for storing the distribution functions (f, feq), velocity components (u, v), and density (rho). Additional arrays include weights (w), lattice velocities (ex, ey), and others for supporting calculations.

Distribution Function Initialization: The distribution functions are initialized to their equilibrium values, considering the initial velocity and density of the fluid.

Collision Process: The core of the LBM, the collision step, is implemented. It involves updating the distribution functions based on local fluid properties, followed by a streaming step where distribution functions are propagated to neighboring nodes. The Bhatnagar-Gross-Krook (BGK) approximation is employed for collision dynamics.

Boundary Conditions:

Bounce-Back: At the cylinder boundary, a bounce-back condition reflects particles, ensuring no flow penetration.
Periodic: Horizontal boundaries are treated with periodic conditions, allowing fluid to re-enter the domain from the opposite side.
Fixed Velocity: At the inlet, a constant velocity is imposed, while the outlet adapts to the flow.
Data Generation for Visualization: The program outputs simulation data in the VTK format, suitable for visualization in tools like Paraview. This feature allows for in-depth analysis of fluid behavior around the cylinder.



Contributing
Contributions to improve the simulation or visualization are welcome. Please feel free to fork the repository, make changes, and submit a pull request.


