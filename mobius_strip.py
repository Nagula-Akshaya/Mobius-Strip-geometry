import numpy as np
from scipy.integrate import simpson
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    def __init__(self, R=1.0, w=0.3, n=200):
        # Initialize Möbius strip parameters:
        # R: radius of the center circle
        # w: width of the strip
        # n: number of mesh points for each parameter (u and v)
        self.R = R
        self.w = w
        self.n = n
        # Create a mesh grid for parameters u (angle) and v (width offset)
        self.u, self.v = np.meshgrid(
            np.linspace(0, 2 * np.pi, n),
            np.linspace(-w / 2, w / 2, n)
        )
        self._generate_mesh()

    def _generate_mesh(self):
        # Parametric equations to generate Möbius strip surface points (x, y, z)
        u, v = self.u, self.v
        self.x = (self.R + v * np.cos(u / 2)) * np.cos(u)
        self.y = (self.R + v * np.cos(u / 2)) * np.sin(u)
        self.z = v * np.sin(u / 2)

    def surface_area(self):
        # Compute partial derivatives with respect to u and v
        dx_du = np.gradient(self.x, self.u[0], axis=1)
        dx_dv = np.gradient(self.x, self.v[:, 0], axis=0)
        dy_du = np.gradient(self.y, self.u[0], axis=1)
        dy_dv = np.gradient(self.y, self.v[:, 0], axis=0)
        dz_du = np.gradient(self.z, self.u[0], axis=1)
        dz_dv = np.gradient(self.z, self.v[:, 0], axis=0)

        # Calculate cross product of the partial derivatives to get surface element normals
        cross_x = dy_du * dz_dv - dz_du * dy_dv
        cross_y = dz_du * dx_dv - dx_du * dz_dv
        cross_z = dx_du * dy_dv - dy_du * dx_dv

        # Compute local surface area element dA from magnitude of normal vector
        dA = np.sqrt(cross_x ** 2 + cross_y ** 2 + cross_z ** 2)

        # Integrate over the surface using 2D Simpson's rule
        area = simpson(simpson(dA, self.v[:, 0], axis=0), self.u[0], axis=0)
        return area

    def edge_length(self):
        # Get coordinates of one boundary edge (v = +w/2)
        x_edge = self.x[-1]
        y_edge = self.y[-1]
        z_edge = self.z[-1]
        total_length = 0

        # Sum Euclidean distances between consecutive points on the edge
        for i in range(1, len(x_edge)):
            p1 = (x_edge[i - 1], y_edge[i - 1], z_edge[i - 1])
            p2 = (x_edge[i], y_edge[i], z_edge[i])
            total_length += euclidean(p1, p2)

        # Close the loop by adding the distance from last to first point
        total_length += euclidean((x_edge[-1], y_edge[-1], z_edge[-1]),
                                  (x_edge[0], y_edge[0], z_edge[0]))

        # Multiply by 2 because the Möbius strip has only one edge that twists back
        return total_length * 2

    def plot(self):
        # Create a 3D plot of the Mobius strip surface
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.x, self.y, self.z, cmap='plasma', edgecolor='k', linewidth=0.1)
        ax.set_title("Möbius Strip")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.tight_layout()
        plt.show()

# Run the script to calculate and visualize the Mobius strip
if __name__ == "__main__":
    mobius = MobiusStrip(R=1.0, w=0.3, n=300)
    print(f"Surface Area ≈ {mobius.surface_area():.4f}")
    print(f"Edge Length ≈ {mobius.edge_length():.4f}")
    mobius.plot()
