import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":

    # ----------------------
    # Create seiscope model
    # ----------------------
    dz = 15
    dx = 15
    nz_new = 301
    nx_new = 501
    zmax = (nz_new - 1) * dz
    xmax = (nx_new - 1) * dx
    extent = [0, xmax, zmax, 0]

    print("Creating Seiscope model...")
    vp = np.zeros(shape=(nz_new, nx_new), dtype=np.float32) + 4.0
    vp[:100, :] = 2.8

    fig = plt.figure(figsize=(12, 3))  # define figure size
    image = plt.imshow(vp, cmap="jet", interpolation='nearest', extent=extent, vmin=1.5, vmax=4.5)
    cbar = plt.colorbar(aspect=10, pad=0.02)
    cbar.set_label('Vp [km/s]', labelpad=10)
    plt.xlabel(r'$x_1$ [m]')
    plt.ylabel(r'$z$ [m]')

    np.savez("Lippmann-Schwinger/Data/seiscope_new_vz_2d.npz", vp)
    savefig_fname = "Lippmann-Schwinger/Fig/seiscope_new_vz_2d.pdf"
    plt.savefig(savefig_fname, format="pdf", bbox_inches="tight", pad_inches=0.01)
    plt.show()

    # --------------------------------------------------------
    # Create perturbation
    # --------------------------------------------------------
    vp_diff = vp * 0.0
    for i in range(100, nz_new):
        for j in range(nx_new):
            if (((i - 100) ** 2 + (j - 250) ** 2) ** 0.5) < 100:
                vp_diff[i, j] = -4.0 + 2.8

    vp_total = vp_diff + vp

    fig = plt.figure(figsize=(12, 3))  # define figure size
    image = plt.imshow(vp_total, cmap="jet", interpolation='nearest', extent=extent, vmin=1.5, vmax=4.5)
    cbar = plt.colorbar(aspect=10, pad=0.02)
    cbar.set_label('Vp [km/s]', labelpad=10)
    plt.xlabel(r'$x_1$ [m]')
    plt.ylabel(r'$z$ [m]')

    np.savez("Lippmann-Schwinger/Data/seiscope_new_2d.npz", vp_total)
    savefig_fname = "Lippmann-Schwinger/Fig/seiscope_new_2d.pdf"
    plt.savefig(savefig_fname, format="pdf", bbox_inches="tight", pad_inches=0.01)
    plt.show()
