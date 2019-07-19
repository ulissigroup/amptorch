import sys
from ase import io
from ase.visualize import view


f = io.read("../datasets/COCu/COCu.traj", ":")
# f = io.read("MD_results/C_Cu_distance.traj", ":")
# f = io.read("MD_results/MLMD_COCu.traj", ":")
f = io.read("MD_results/MLMD_LJ_COCu.traj", ":")
images = []
for i in range(100):
    images.append(f[i])

gif = []
for i, s in enumerate(images):
    for k, v in enumerate(s):
        s[k].position[2] -= 10
    s.rotate(-90, "x", rotate_cell=True)
    s.rotate(45, "y", rotate_cell=True)
    s.rotate(20, "x", rotate_cell=True)
    # s.rotate(45, "y", rotate_cell=True)
    # s.rotate(20, "x", rotate_cell=True)
    # s.set_cell((1, 1, 1))
    s = s.repeat((3, 3, 1))
    s.set_cell((0, 0, 0))
    gif.append(s)

for idx, image in enumerate(gif):
    name = str(idx)+'.pov'
    io.write(name, image)
# io.write("MD_results/movies/COCu_ML.gif", gif, interval=100)
