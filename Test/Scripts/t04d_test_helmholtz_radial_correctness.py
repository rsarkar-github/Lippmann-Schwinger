import sys
import numpy as np
import time
from scipy.sparse.linalg import splu, lsqr
from matplotlib import pyplot as plt
from ...Solver.HelmholtzOperators import create_helmholtz3d_matrix
from ...Solver.ScatteringIntegralConstantVel import TruncatedKernelConstantVel3d


if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) < 2:
        raise ValueError("Program missing command line arguments.")
    plot_flag = bool(int(sys.argv[1]))

    # Lippmann-Schwinger solver
    n = 201
    d = 1.0 / (n - 1)
    precision = np.complex64
    f = 15.0
    v0 = 1.0
    k = 2 * np.pi * f / v0

    op = TruncatedKernelConstantVel3d(n=n, k=k, precision=precision)

    u = np.zeros(shape=(n, n, n), dtype=precision)
    u[int(n/2), int(n/2), int(n/2)] = 1.0
    sol1 = u * 0

    start_t = time.time()
    op.convolve_kernel(u=u, output=sol1)
    end_t = time.time()
    print("Total time to execute convolution: ", "{:4.2f}".format(end_t - start_t), " s \n")

    # Helmholtz solver
    pml_cells = 10
    n1 = 221
    vel = np.zeros(shape=(n1, n1, n1), dtype=np.float32) + v0
    omega = 2 * np.pi * f
    a = d * (n1 - 1)

    mat_3d = create_helmholtz3d_matrix(
        a1=a,
        a2=a,
        a3=a,
        pad1=pml_cells,
        pad2=pml_cells,
        pad3=pml_cells,
        omega=omega,
        precision=precision,
        vel=vel,
        pml_damping=50.0,
        adj=False,
        warnings=False
    )
    u = np.zeros(shape=(n1, n1, n1), dtype=precision)
    u[int(n1 / 2), int(n1 / 2), int(n1 / 2)] = 1.0

    tol = 1e-3
    sol, istop, itn, r1norm = lsqr(
        mat_3d,
        np.reshape(u, newshape=(n1 * n1 * n1, 1)),
        atol=0,
        btol=tol,
        show=True,
        iter_lim=100
    )[:4]
    print(itn, r1norm)
    sol2 = np.reshape(sol, newshape=(n1, n1, n1))[10:211, 10:211, 10:211]

    print("Relative error = ", np.linalg.norm(sol1 - sol2) / np.linalg.norm(sol1))

    scale = 1e-4
    plt.imshow(np.real(sol1[:, :, int(n/2)]), cmap="Greys", vmin=-scale, vmax=scale)
    plt.show()

    plt.imshow(np.real(sol2[:, :, int(n / 2)]), cmap="Greys", vmin=-scale, vmax=scale)
    plt.show()
