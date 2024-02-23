import time
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.managers import SharedMemoryManager
from scipy.sparse.linalg import LinearOperator
from Solver.ScatteringIntegralGeneralVz import TruncatedKernelGeneralVz2d


if __name__ == "__main__":

    # Load files
    filepath = "Lippmann-Schwinger/Data/p04c-seiscope-new-vz-2d.npz"
    filepath1 = "Lippmann-Schwinger/Data/p04c-seiscope-new-2d.npz"
    filepath2 = "Lippmann-Schwinger/Data/p06-seiscope-source.npz"

    # Freq and solver
    freq = 5.0   # in Hz
    solver_name = "gmres"

    # ----------------------------------------------
    # Load vz and calculate psi
    # Load initial solution
    # ----------------------------------------------
    with np.load(filepath) as data:
        vel = data["arr_0"]
    vel_trace = vel[:, 0]
    n1_vel_trace = vel_trace.shape[0]
    vel_trace = np.reshape(vel_trace, newshape=(n1_vel_trace, 1)).astype(np.float32)

    with np.load(filepath1) as data:
        vel1 = data["arr_0"]

    psi_ = (1.0 / vel) ** 2.0 - (1.0 / vel1) ** 2.0

    # ----------------------------------------------
    # Set parameters
    # ----------------------------------------------
    n_ = 351
    nz_ = 251
    a_ = 0.
    b_ = a_ + (1.0 / (n_ - 1)) * (nz_ - 1)
    freq_ = freq * 5.25
    omega_ = 2 * np.pi * freq_
    m_ = 4
    sigma_ = 0.0015
    precision_ = np.complex64
    green_func_dir_ = filepath = "Lippmann-Schwinger/Data/p05-green-func-2-0"
    num_threads_ = 4
    vz_ = np.zeros(shape=(nz_, 1), dtype=np.float32) + vel_trace

    psi_ = psi_.astype(precision_)

    # ----------------------------------------------
    # Load source
    # ----------------------------------------------
    with np.load(filepath2) as data:
        sou_ = data["arr_0"]

    sou_ = sou_.astype(precision_)

    # ----------------------------------------------
    # Get Green's function
    # ----------------------------------------------
    op = TruncatedKernelGeneralVz2d(
        n=n_,
        nz=nz_,
        a=a_,
        b=b_,
        k=omega_,
        vz=vz_,
        m=m_,
        sigma=sigma_,
        precision=precision_,
        green_func_dir=green_func_dir_,
        num_threads=num_threads_,
        verbose=False,
        light_mode=True
    )
    op.set_parameters(
        n=n_,
        nz=nz_,
        a=a_,
        b=b_,
        k=omega_,
        vz=vz_,
        m=m_,
        sigma=sigma_,
        precision=precision_,
        green_func_dir=green_func_dir_,
        num_threads=num_threads_,
        verbose=False
    )

    green_func = op.greens_func

    # Calculate number of bytes
    num_bytes = self._nz * self._nz * self._num_bins_non_neg * self._num_bins_non_neg
    if self._precision == np.complex64:
        num_bytes *= 8
    if self._precision == np.complex128:
        num_bytes *= 16

    # ----------------------------------------------
    # Setup multiprocessing workflow
    # ----------------------------------------------

    # Read in multiprocessing mode
    with SharedMemoryManager() as smm:


        # Create shared memory
        sm = smm.SharedMemory(size=num_bytes)
        self._green_func = ndarray(
            shape=(self._nz, self._nz, self._num_bins_non_neg, self._num_bins_non_neg),
            dtype=self._precision,
            buffer=sm.buf
        )
        self._green_func *= 0

        param_tuple_list = [
            (
                self._green_func_dir,
                ii,
                self._nz,
                self._n,
                self._precision,
                sm.name
            ) for ii in range(self._nz)
        ]

        print("\nReading Green's function...")

        with Pool(min(len(param_tuple_list), mp.cpu_count(), self._num_threads)) as pool:
            max_ = len(param_tuple_list)

            with tqdm(total=max_) as pbar:
                for _ in pool.imap_unordered(func_read3D, param_tuple_list):
                    pbar.update()


    rhs_ = np.zeros(shape=(nz_, n_), dtype=precision_)
    t1 = time.time()
    op.apply_kernel(u=sou_, output=rhs_)
    t2 = time.time()
