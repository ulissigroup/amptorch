import sys
import os
from ase import io
from ase.visualize import view
from ase.build import add_adsorbate, molecule
from ase import Atoms


def create_movie(location, destination, filename):
    if not os.path.exists(destination):
        os.mkdir(destination)

    data = io.read(location + filename + ".traj", ":")
    images = [data[i] for i in range(100)]

    gif = []
    for i, s in enumerate(images):
        # for k, v in enumerate(s):
        # s[k].position[2] -= 10
        s.rotate(-45, "z", rotate_cell=True)
        # del s[[atom.index for atom in s if atom.position[1]>=15]]
        s.rotate(-70, "x", rotate_cell=True)
        s = s.repeat((3, 3, 1))
        # s.extend(Atoms("CO", [(-17.690, 17, 40), (-17.690, 18.1, 40)]))
        s.wrap()
        # s.wrap(center=(1, 1, 1))
        s.set_cell((0, 0, 0))
        gif.append(s)

    io.write(destination + filename + ".gif", gif, interval=100)


def retrieve_files(prefix, num_images):
    files = []
    for i in range(num_images):
        file_ML = "".join([prefix, "_%s" % (i + 1)])
        file_LJ = "".join([prefix, "_LJ", "_%s" % (i + 1)])
        files.append(file_ML)
        files.append(file_LJ)
    return files


location = "MD_results/COCu/pbc_300K/l2amp/paper/"
# location = "../../datasets/COCu/"
destination = "MD_results/movies/COCu/pbc_300K/l2amp/"
# files = retrieve_files("MLMD_COCu_pbc_300K_log", 3)

files = [
    "MLMD_COCu_pbc_300K_l2amp_8SF_1",
    "MLMD_COCu_pbc_300K_l2amp_8SF_LJ_1",
    "MLMD_COCu_pbc_300K_l2amp_8SF_5_resample_1",
    "MLMD_COCu_pbc_300K_l2amp_8SF_5_resample_2",
    "MLMD_COCu_pbc_300K_l2amp_8SF_LJ_5_resample_1",
    "MLMD_COCu_pbc_300K_l2amp_8SF_LJ_5_resample_2",
]

for file in files:
    create_movie(location, destination, file)
