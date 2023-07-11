
import polyscope as ps
import numpy as np
import openmesh as om
import argparse
import os
from implementation import *

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--file", type=str, default="Meshes/Geometry/cube.off", help="Nombre del archivo")
  filename = parser.parse_args().file

  mesh = om.read_polymesh(filename)

  ps.init()
  ps_mesh_original = ps.register_surface_mesh(
    "mesh_original", 
    mesh.points(), 
    mesh.face_vertex_indices(),
    transparency=0.5
    )
  
  ps_point_cloud = ps.register_point_cloud(
    "point_cloud", 
    catmull_clark(mesh).points(),
    radius=0.01
    )
  ps.show()

if __name__ == "__main__": main()