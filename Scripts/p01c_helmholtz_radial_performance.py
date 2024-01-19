import sys
import numpy as np
import time
from scipy import interpolate
from scipy.sparse.linalg import gmres
from matplotlib import pyplot as plt
from ..Solver.HelmholtzOperators import create_helmholtz2d_matrix_radial
from ..Utilities import TypeChecker
from ..Utilities.Utils import make_velocity_from_trace, make_velocity3d_from_trace, extend_vel_trace_1d



if __name__ == "__main__":

    # Check arguments
    if len(sys.argv) < 2:
        raise ValueError("Program missing command line arguments.")

    freq = float(sys.argv[1])

    # Load Marmousi velocity trace
    with np.load("Lippmann-Schwinger/Data/marmousi-vp-vz.npz") as data:
        vel_trace = data["arr_0"]
    vel_trace /= 1000.0
    n1_vel_trace = vel_trace.shape[0]
    vel_trace = np.reshape(vel_trace, newshape=(n1_vel_trace, 1))

    # Define frequency, calculate min & max wavelength
    omega = freq * 2 * np.pi
    precision = np.complex128

    vmin = np.min(vel_trace)
    vmax = np.max(vel_trace)
    lambda_min = vmin / freq
    lambda_max = vmax / freq

    # Set grid extent, & calculate minimum grid spacing
    delta_base = 0.015
    a1 = n1_vel_trace * delta_base  # in km (should not change)
    a2 = 400 * delta_base           # in km (should not change)

    # Grid refining to do (different experiments)
    fac = [1, 2, 4]

    # Function for creating all objects for the experiment
    def create_objects_for_experiment(factor):
        """
        :param factor: factor by which to refine the computational grid
        :return: New grid parameters, velocity
        """

        # New grid spacing
        delta = delta_base / factor

        # Get new extents
        a1_pad, a2_pad, pad1_cells, pad2_cells = make_grid_params(
            a1_=a1,
            a2_=a2,
            delta_=delta,
            lambda_min_=lambda_min,
            lambda_max_=lambda_max
        )

        # Calculate grid points in each direction
        n1_ = int(np.round(a1_pad / delta) + 1)
        n2_ = int(np.round(a2_pad / delta) + 1)
        n1_vel = int(np.round(a1 / delta) + 1)

        # Modify vel_trace
        vel_trace_mod = make_velocity_from_trace(vel_trace_=vel_trace, n1_=n1_vel, n2_=1)
        vel_trace_mod = extend_vel_trace_1d(vel_trace_=vel_trace_mod, pad_cells_=pad1_cells)
        vel = make_velocity_from_trace(vel_trace_=vel_trace_mod, n1_=n1_, n2_=n2_)

        # Create source centered at original grid (x2 = 0, x1 = a1 / 10, gaussian std = delta)
        x1 = np.linspace(start=0, stop=a1_pad, num=n1_, endpoint=True)
        x2 = np.linspace(start=0, stop=a2_pad, num=n2_, endpoint=True)
        x1v, x2v = np.meshgrid(x1, x2)
        sigma = delta_base
        source = np.exp((-1) * ((x1v - a1_pad / 10) ** 2 + x2v ** 2) / (2 * sigma * sigma))
        source = source.T

        return a1_pad, a2_pad, pad1_cells, pad2_cells, vel, source

    for i, item in enumerate(fac):

        print("\n\n---------------------------------------------------------")

        a1_full, a2_full, pad1, pad2, vel_array, src = create_objects_for_experiment(factor=item)

        mat = create_helmholtz2d_matrix_radial(
            a1=a1_full,
            a2=a2_full,
            pad1=pad1,
            pad2=pad2,
            omega=omega,
            precision=precision,
            vel=vel_array,
            pml_damping=100.0,
            adj=False,
            warnings=True
        )

        n1, n2 = vel_array.shape
        print("n1 = ", n1, "n2 = ", n2)

        def make_callback():
            closure_variables = dict(counter=0, residuals=[])

            def callback(residuals):
                closure_variables["counter"] += 1
                closure_variables["residuals"].append(residuals)
                print(closure_variables["counter"], residuals)

            return callback

        # GMRES
        tol = 1e-6
        t_start = time.time()
        sol, exitcode = gmres(
            mat,
            np.reshape(src, newshape=(n1 * n2, 1)),
            maxiter=50000,
            restart=200,
            callback=make_callback()
        )
        print("\nExitcode ", exitcode)
        t_end = time.time()
        print("Time to solve GMRES (restart = 200) is ", "{:4.2f}".format(t_end - t_start), "s")
        sol = np.reshape(sol, newshape=(n1, n2))

        scale = 1e-5
        fig, ax = plt.subplots(1, 1)
        im = ax.imshow(np.real(sol), cmap="Greys", vmin=-scale, vmax=scale)
        plt.show()
