import numpy as np
import time
from scipy.sparse.linalg import gmres
from matplotlib import pyplot as plt
from ...Solver.HelmholtzOperators import create_helmholtz3d_matrix
from ...Solver.ScatteringIntegralGeneralVz import TruncatedKernelGeneralVz3d



if __name__ == "__main__":

    sigma = 0.03
    n = 201
    d = 1.0 / (n - 1)
    precision = np.complex64
    f = 10.0
    v0 = 1.0
    omega = 2 * np.pi * f

    # Compute Lippmann-Schwinger solver
    a = 0.
    b = a + (1.0 / (n - 1)) * (n - 1)
    m = 5
    vz = np.zeros(shape=(n, 1), dtype=np.float32) + v0

    def lippmann_schwinger():
        op = TruncatedKernelGeneralVz3d(
            n=n,
            nz=n,
            a=a,
            b=b,
            k=omega,
            vz=vz,
            m=m,
            sigma=3 * d/m,
            precision=precision,
            green_func_dir="Lippmann-Schwinger/Test/Data/t13",
            num_threads=50,
            verbose=False,
            light_mode=False
        )

        grid = np.linspace(start=-0.5, stop=0.5, num=n, endpoint=True)
        x1, x2, z1 = np.meshgrid(grid, grid, grid, indexing="ij")
        u = np.exp(-1.0 * (x1 ** 2 + x2 ** 2 + z1 ** 2) / (2 * sigma * sigma))
        u = u.astype(precision)
        sol1 = u * 0

        start_t = time.time()
        op.apply_kernel(u=u, output=sol1)
        end_t = time.time()
        print("Total time to execute convolution: ", "{:4.2f}".format(end_t - start_t), " s \n")

        return sol1, u

    sol1, u = lippmann_schwinger()

    def make_callback():
        closure_variables = dict(counter=0, residuals=[])

        def callback(residuals):
            closure_variables["counter"] += 1
            closure_variables["residuals"].append(residuals)
            print(closure_variables["counter"], residuals)

        return callback

    def helmholtz():

        pml_cells = 10
        n1 = n + 2 * pml_cells
        vel = np.zeros(shape=(n1, n1, n1), dtype=np.float32) + v0
        a1 = d * (n1 - 1)

        mat = create_helmholtz3d_matrix(
            a1=a1,
            a2=a1,
            a3=a1,
            pad1=pml_cells,
            pad2=pml_cells,
            pad3=pml_cells,
            omega=omega,
            precision=precision,
            vel=vel,
            pml_damping=50.0,
            adj=False,
            warnings=True
        )
        u1 = np.zeros(shape=(n1, n1, n1), dtype=precision)
        u1[pml_cells:n+pml_cells, pml_cells:n+pml_cells, pml_cells:n+pml_cells] = u

        tol = 1e-5

        # GMRES
        sol, exitcode = gmres(
            mat,
            np.reshape(u1, newshape=(n1 * n1 * n1, 1)),
            maxiter=10000,
            restart=20,
            callback=make_callback(),
            tol=tol
        )
        print("\nExitcode ", exitcode)

        sol = np.reshape(sol, newshape=(n1, n1, n1))
        sol2 = np.reshape(sol, newshape=(n1, n1, n1))[
               pml_cells:n+pml_cells,
               pml_cells:n+pml_cells,
               pml_cells:n+pml_cells
        ]

        return sol2

    sol2 = helmholtz()

    print("Relative error = ", np.linalg.norm(sol1 - sol2) / np.linalg.norm(sol2))
