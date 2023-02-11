using HDF5
using Interpolations

DATA_DIR = "data/"
DATA_PATH = "jupiter/outGridVelocity.h5"

function load_data()
    fid = h5open(DATA_DIR * DATA_PATH, "r")

    u = read(fid, "vx")
    v = read(fid, "vy")
    return 
end