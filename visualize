import numpy as np
import matplotlib.pyplot as plt

class LatticeBoltzmann:
    def __init__(self, nx=400, ny=100, tau=0.6, u_lb=0.04, n_steps=1000):
        self.nx = nx
        self.ny = ny
        self.tau = tau
        self.u_lb = u_lb
        self.n_steps = n_steps
        self.q = 9
        self.w = [4/9] + [1/9]*4 + [1/36]*4
        self.cx = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
        self.cy = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
        self.opp = [0, 3, 4, 1, 2, 7, 8, 5, 6]

    def initialize(self):
        rho = np.ones((self.nx, self.ny))
        ux = np.zeros((self.nx, self.ny))
        uy = np.zeros((self.nx, self.ny))
        feq = np.zeros((self.q, self.nx, self.ny))
        for i in range(self.q):
            cu = 3 * (self.cx[i]*ux + self.cy[i]*uy)
            feq[i, :, :] = rho * self.w[i] * (1 + cu + 0.5*cu**2 - 1.5*(ux**2 + uy**2))
        return feq

    def collision_and_streaming(self, f):
        rho = np.sum(f, axis=0)
        ux = np.sum(f * self.cx[:, None, None], axis=0) / rho
        uy = np.sum(f * self.cy[:, None, None], axis=0) / rho

        for i in range(self.q):
            cu = 3 * (self.cx[i]*ux + self.cy[i]*uy)
            feq = rho * self.w[i] * (1 + cu + 0.5*cu**2 - 1.5*(ux**2 + uy**2))
            f[i, :, :] += -(1.0/self.tau) * (f[i, :, :] - feq)

        for i in range(self.q):
            f[i, :, :] = np.roll(f[i, :, :], self.cx[i], axis=0)
            f[i, :, :] = np.roll(f[i, :, :], self.cy[i], axis=1)

        return f, ux, uy

    def apply_boundary_conditions(self, f):
        cylinder = (self.nx//4, self.ny//2, self.ny//10)
        mask = (np.arange(self.nx)[:, None] - cylinder[0])**2 + (np.arange(self.ny)[None, :] - cylinder[1])**2 <= cylinder[2]**2

        for i in range(self.q):
            f[self.opp[i], mask] = f[i, mask]

        f[:, 0, :] = self.initialize()[:, 0, :]  # Inlet
        f[:, -1, :] = f[:, -2, :]  # Outlet
        return f

    def plot(self, ux, uy, time_step):
        plt.imshow(np.sqrt(ux**2 + uy**2).transpose(), origin='lower', cmap='viridis')
        plt.colorbar()
        plt.title(f"Velocity field at timestep {time_step}")
        plt.show()


def main():
    lbm = LatticeBoltzmann()
    f = lbm.initialize()

    for time_step in range(lbm.n_steps):
        f, ux, uy = lbm.collision_and_streaming(f)
        f = lbm.apply_boundary_conditions(f)

        if time_step % 100 == 0:
            lbm.plot(ux, uy, time_step)


if __name__ == "__main__":
    main()
