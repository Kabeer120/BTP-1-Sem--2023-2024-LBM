import numpy as np
import matplotlib.pyplot as plt

# Define the LatticeBoltzmann class
class LatticeBoltzmann:
    # Initialize the class with default parameters
    def __init__(self, nx=400, ny=100, tau=0.6, u_lb=0.04, n_steps=1000):
        self.nx = nx  # Number of lattice nodes in the x-direction
        self.ny = ny  # Number of lattice nodes in the y-direction
        self.tau = tau  # Relaxation time
        self.u_lb = u_lb  # Lattice unit for velocity
        self.n_steps = n_steps  # Number of simulation time steps
        self.q = 9  # Number of discrete velocity directions in the D2Q9 model
        self.w = [4/9] + [1/9]*4 + [1/36]*4  # D2Q9 weights for each direction
        self.cx = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])  # x components of velocity directions
        self.cy = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])  # y components of velocity directions
        self.opp = [0, 3, 4, 1, 2, 7, 8, 5, 6]  # Opposite direction indices for bounce-back

    # Initialize density and velocity fields
    def initialize(self):
        rho = np.ones((self.nx, self.ny))  # Initialize density to 1 for all lattice points
        ux = np.zeros((self.nx, self.ny))  # Initialize x-velocity to 0 for all lattice points
        uy = np.zeros((self.nx, self.ny))  # Initialize y-velocity to 0 for all lattice points
        feq = np.zeros((self.q, self.nx, self.ny))  # Initialize distribution functions
        for i in range(self.q):  # Calculate equilibrium distribution function for each direction
            cu = 3 * (self.cx[i]*ux + self.cy[i]*uy)
            feq[i, :, :] = rho * self.w[i] * (1 + cu + 0.5*cu**2 - 1.5*(ux**2 + uy**2))
        return feq

    # Perform collision and streaming steps
    def collision_and_streaming(self, f):
        rho = np.sum(f, axis=0)  # Compute density as the sum of distribution functions
        ux = np.sum(f * self.cx[:, None, None], axis=0) / rho  # Compute x-velocity
        uy = np.sum(f * self.cy[:, None, None], axis=0) / rho  # Compute y-velocity

        for i in range(self.q):  # Collision step: calculate and apply changes to distribution functions
            cu = 3 * (self.cx[i]*ux + self.cy[i]*uy)
            feq = rho * self.w[i] * (1 + cu + 0.5*cu**2 - 1.5*(ux**2 + uy**2))
            f[i, :, :] += -(1.0/self.tau) * (f[i, :, :] - feq)

        for i in range(self.q):  # Streaming step: shift distribution functions along their velocity directions
            f[i, :, :] = np.roll(f[i, :, :], self.cx[i], axis=0)
            f[i, :, :] = np.roll(f[i, :, :], self.cy[i], axis=1)

        return f, ux, uy

    # Apply boundary conditions to the lattice
    def apply_boundary_conditions(self, f):
        # Define the cylinder boundary in the lattice
        cylinder = (self.nx//4, self.ny//2, self.ny//10)
        # Create a mask for the cylinder
        mask = (np.arange(self.nx)[:, None] - cylinder[0])**2 + (np.arange(self.ny)[None, :] - cylinder[1])**2 <= cylinder[2]**2

        for i in range(self.q):  # Apply bounce-back boundary condition
            f[self.opp[i], mask] = f[i, mask]

        # Set inlet and outlet boundary conditions
        f[:, 0, :] = self.initialize()[:, 0, :]  # Inlet condition
        f[:, -1, :] = f[:, -2, :]  # Outlet condition
        return f

    # Generate and display the plot of the velocity field
    def plot(self, ux, uy, time_step):
        plt.imshow(np.sqrt(ux**2 + uy**2).transpose(), origin='lower', cmap='viridis')  # Plot the velocity magnitude
        plt.colorbar()  # Add a color bar to indicate the magnitude
        plt.title(f"Velocity field at timestep {time_step}")  # Set the title of the plot
        plt.show()  # Display the plot

# Main function to execute the simulation
def main():
    lbm = LatticeBoltzmann()  # Create an instance of the LatticeBoltzmann class
    f = lbm.initialize()  # Initialize the simulation

    for time_step in range(lbm.n_steps):  # Main simulation loop
        f, ux, uy = lbm.collision_and_streaming(f)  # Perform collision and streaming
        f = lbm.apply_boundary_conditions(f)  # Apply boundary conditions

        if time_step % 100 == 0:  # Plot every 100 time steps
            lbm.plot(ux, uy, time_step)

# Entry point of the script
if __name__ == "__main__":
    main()  # Call the main function
