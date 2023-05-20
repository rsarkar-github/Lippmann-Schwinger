import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
import time
from ...Solver.HelmholtzOperators import create_helmholtz2d_matrix


if __name__ == "__main__":

    a1_ = 1.0
    a2_ = 1.0

    freq = 15.0
    omega = freq * 2 * np.pi
    precision = np.complex64

    v0 = 2.0
    lambda_ = v0 / freq
    grid_spacing_min = lambda_ / 10.0
    padwidth_min = lambda_

    n1_start = int(a1_ / grid_spacing_min) + 1
    n2_start = int(a2_ / grid_spacing_min) + 1

    def make_odd(ii):
        if ii % 2 == 1:
            return ii
        else:
            return ii + 1

    n1_start = make_odd(n1_start)
    n2_start = make_odd(n2_start)

    fac = [1.0, 2.0, 4.0, 6.0]
    arr_n1 = [int(n1_start * item) for item in fac]
    arr_n2 = [int(n2_start * item) for item in fac]

    a1_arr = []
    a2_arr = []
    pad1_arr = []
    pad2_arr = []
    n1_total_arr = []
    n2_total_arr = []

    for i, item in enumerate(fac):

        dx1 = a1_ / arr_n1[i]
        dx2 = a2_ / arr_n2[i]

        pad1_arr.append(int(padwidth_min / dx1) + 1)
        pad2_arr.append(int(padwidth_min / dx2) + 1)
        n1_total_arr.append(arr_n1[i] + 2 * pad1_arr[i])
        n2_total_arr.append(arr_n2[i] + 2 * pad2_arr[i])
        a1_arr.append(dx1 * n1_total_arr[i])
        a2_arr.append(dx2 * n2_total_arr[i])


    def factorize(n1, n2, pad1, pad2, a1, a2):

        vel = np.zeros(shape=(n1, n2), dtype=np.float32)
        vel += v0

        # Create Helmholtz matrix
        mat = create_helmholtz2d_matrix(
            a1=a1,
            a2=a2,
            pad1=pad1,
            pad2=pad2,
            omega=omega,
            precision=precision,
            vel=vel,
            pml_damping=50.0,
            adj=False,
            warnings=True
        )

        start_t = time.time()
        mat_lu = splu(mat)
        end_t = time.time()
        print("\nTotal time to LU factorize: ", "{:4.2f}".format(end_t - start_t), " s \n")

        return end_t - start_t

    fac_times = []
    for i, item in enumerate(fac):

        t = factorize(
            n1=n1_total_arr[i],
            n2=n2_total_arr[i],
            pad1=pad1_arr[i],
            pad2=pad2_arr[i],
            a1=a1_arr[i],
            a2=a2_arr[i]
        )

        fac_times.append(t)
        print("n1 = ", n1_total_arr[i], ", n2 = ", n2_total_arr[i], ", time = ", "{:4.2f}".format(t))
