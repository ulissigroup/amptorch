import sys
from ase import io
from ase.visualize import view
from ase.build import add_adsorbate, molecule
from ase import Atoms


# data = io.read("../datasets/COCu/COCu_pbc.traj", ":")
data = io.read("../datasets/COCu/COCu.traj", ":")
# data  = io.read("MD_results/MLMD_COCu.traj", ":")
# data = io.read("MD_results/MLMD_COCu_pbc.traj", ":")
# data = io.read("MD_results/MLMD_LJ_COCu_pbc.traj", ":")
# data = io.read("MD_results/MLMD_LJ_COCu.traj", ":")
# f = io.read("../datasets/surface_water.traj", ":")
# data = io.read("MD_results/C_Cu_distance.traj", ":")
images = []
for i in range(100):
    images.append(data[i])

# for i, image in enumerate(images):
    # image = image.repeat((2, 2, 2))
    # image.wrap()
    # image.write('images/%04d.pov'%i, run_povray=True, transparent=False,
            # pause=False, canvas_width=400, rotation='-45z, -70x',
            # show_unit_cell=2)
# !ffmpeg -i images/%04d.png output.gif -y

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
    s.set_cell((0,0,0))
    gif.append(s)

# for idx, image in enumerate(gif):
    # name = str(idx)+'.pov'
    # io.write(name, image)
# io.write("MD_results/movies/COCu_ML_wrap.gif", gif, interval=100)
# io.write("test.gif", f, interval=100)
