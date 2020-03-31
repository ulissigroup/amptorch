import sys
import os
from ase import io
from ase.visualize import view
from ase.build import add_adsorbate, molecule
from ase import Atoms


def create_movie(location, destination, filename):
    os.makedirs(destination, exist_ok=True)

    data = io.read(location + filename + ".traj", ":500")

    gif = []
    for i, s in enumerate(data):
        s.rotate(-45, "z", rotate_cell=True)
        s.rotate(-70, "x", rotate_cell=True)
        s = s.repeat((3, 3, 1))
        s.wrap()
        s.set_cell((0, 0, 0))
        gif.append(s)

    io.write(destination + filename + ".gif", gif, interval=10)


def retrieve_files(prefix, num_images):
    files = []
    for i in range(num_images):
        file_ML = "".join([prefix, "_%s" % (i + 1)])
        file_LJ = "".join([prefix, "_LJ", "_%s" % (i + 1)])
        files.append(file_ML)
        files.append(file_LJ)
    return files


# location = "../ber_results/paper/"
location = "../../datasets/COCu/"
destination = "./"
# destination = "../ber_results/paper/animations/"
# files = retrieve_files("MLMD_COCu_pbc_300K_log", 3)

files = [
        "COCu_aimd_300K"
        # "COCu_ber_2ps_optim_fixmd_LJ_100_iter_1",
        # "COCu_ber_2ps_optim_fixmd_LJ_100_iter_2",
        # "COCu_ber_2ps_optim_fixmd_LJ_100_iter_3",
        # "COCu_ber_2ps_optim_fixmd_LJ_100_iter_4",
        # "COCu_ber_2ps_optim_fixmd_LJ_100_iter_5"
]

for file in files:
    create_movie(location, destination, file)
    print("Movie created: %s" % file)


