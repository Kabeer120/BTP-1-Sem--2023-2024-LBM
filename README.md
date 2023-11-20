# BTP-1-Sem--2023-2024-LBM
This Repository contains the code files for Lattice Boltzmann Method for our BTP Project.

# Lattice Boltzmann Method for Flow Past a Cylinder

## Introduction
This project simulates the fluid flow past a cylinder using the Lattice Boltzmann Method (LBM), a computational fluid dynamics technique. It consists of two main components: a C++ simulation file (`flowPastCylinder.cpp`) and a Python visualization script (`visualize.py`).

## Features
- 2D fluid flow simulation around a cylindrical obstacle.
- Visualization of velocity fields and vorticity.
- Customizable parameters for different flow scenarios.

## Dependencies
- C++ Compiler (e.g., g++, clang)
- Python 3.x
- NumPy
- Matplotlib

## Getting Started

To use this simulation, first clone the repository to your local machine:

```bash
git clone https://github.com/Kabeer120/BTP-1-Sem--2023-2024-LBM.git
cd BTP-1-Sem--2023-2024-LBM
'''
Ensure you have a g++ compiler

g++ -o flowPastCylinder flowPastCylinder.cpp -std=c++11

Visualize from  -

python visualize.py
