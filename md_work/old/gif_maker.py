import os
from os import listdir
from ase import io
from ase import Atoms


def create_movie(location, destination, filename, saveas):
    if not os.path.exists(destination):
        os.mkdir(destination)

    data = io.read(location + filename, ":")
    io.write(destination + saveas + ".gif", data, interval=2)


# dir of traj files to be made to gifs
location = "../../datasets/ktran/"
# where to save gifs
destination = "./"

# files = listdir(location)
files = ["OUTCAR"]

for file in files:
    create_movie(location, destination, file, saveas='test')
    print("Movie created: %s" % file)
