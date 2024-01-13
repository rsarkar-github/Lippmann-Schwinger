import numpy as np
import time
from scipy.sparse.linalg import splu
from scipy import special
from matplotlib import pyplot as plt
from ...Solver.HelmholtzOperators import create_helmholtz2d_matrix_radial
from ...Solver.ScatteringIntegralConstantVelStorageOptimized import TruncatedKernelConstantVel3d


# TODO: problems detected with Helmholtz radial solvers

if __name__ == "__main__":

    # Lippmann-Schwinger solver
    n = 201
    d = 1.0 / (n - 1)
    precision = np.complex64
    f = 10.0
    v0 = 1.0
    k = 2 * np.pi * f / v0
    sigma = 0.05

    def create_source():

        grid = np.linspace(start=-0.5, stop=0.5, num=n, endpoint=True)
        x1, y1, z1 = np.meshgrid(grid, grid, grid, indexing="ij")
        u = np.exp(-1.0 * (x1 ** 2 + y1 ** 2 + z1 ** 2) / (2 * sigma * sigma))
        u *= (1.0 / (sigma ** 3)) / ((2 * np.pi) ** 1.5)
        u = u.astype(precision)

        return u

    u = create_source()

    def helmholtz():

        pml_cells = 10
        n1 = 221
        n2 = 111
        vel = np.zeros(shape=(n1, n2), dtype=np.float32) + v0
        omega = 2 * np.pi * f
        a1 = d * (n1 - 1)
        a2 = d * (n2 - 1)

        mat = create_helmholtz2d_matrix_radial(
            a1=a1,
            a2=a2,
            pad1=pml_cells,
            pad2=pml_cells,
            omega=omega,
            precision=precision,
            vel=vel,
            pml_damping=50.0,
            adj=False,
            warnings=True
        )
        u1 = np.zeros(shape=(n1, n2), dtype=precision)
        u1[10:211, 0:101] = u[:, 100:, int(n/2)]

        # plt.imshow(np.real(u1), cmap="Greys")
        plt.imshow(np.real(u[:, :, int(n/2)]), cmap="Greys")
        plt.show()

        mat_lu = splu(mat)
        sol = mat_lu.solve(np.reshape(u1, newshape=(n1 * n2, 1)))
        sol2 = np.reshape(sol, newshape=(n1, n2))[10:211, 0:101]

        return sol2

    sol2 = helmholtz()
    scale = 1e-3
    plt.imshow(np.real(sol2), cmap="Greys", vmin=-scale, vmax=scale)
    plt.show()

    def analytical():

        sol_true = np.zeros(shape=(n, n, n), dtype=precision)

        j = complex(0, 1)
        for i1 in range(n):
            for i2 in range(n):
                for i3 in range(n):

                    coord_x = -0.5 + i1 * d
                    coord_y = -0.5 + i2 * d
                    coord_z = -0.5 + i3 * d
                    r = (coord_x ** 2 + coord_y ** 2 + coord_z ** 2) ** 0.5

                    if r <= 1e-6:
                        sol_true[i1, i2, i3] = 1.0

                    else:
                        t = ((sigma ** 2) * k * j - r) / (sigma * np.sqrt(2))
                        t = special.erf(t)
                        t = t * np.exp(-j * k * r)
                        t = np.real(t) - j * np.sin(k * r)
                        t = t * np.exp(-0.5 * k * k * sigma * sigma) / (4 * np.pi * r)
                        sol_true[i1, i2, i3] = t

        return sol_true

    sol_true = analytical()
    sol_true = sol_true[:, int(n / 2):, int(n / 2)]

    plt.imshow(np.real(sol_true), cmap="Greys", vmin=-scale, vmax=scale)
    plt.show()

    def lippmann_schwinger():

        op = TruncatedKernelConstantVel3d(n=n, k=k, precision=precision)
        sol = u * 0

        start_t = time.time()
        op.convolve_kernel(u=u, output=sol)
        end_t = time.time()
        print("Total time to execute convolution: ", "{:4.2f}".format(end_t - start_t), " s \n")

        return sol

    sol1 = lippmann_schwinger()
    sol1 = sol1[:, int(n / 2):, int(n/2)]

    plt.imshow(np.real(sol1), cmap="Greys", vmin=-scale, vmax=scale)
    plt.show()

    print("Relative error = ", np.linalg.norm(sol1 - sol_true) / np.linalg.norm(sol1))
