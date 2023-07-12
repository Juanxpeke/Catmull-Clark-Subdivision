
import polyscope as ps
import openmesh as om
import argparse
from implementation import *

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--file", type=str, default="Meshes/Geometry/cube.off", help="Nombre del archivo")
  parser.add_argument("--iter", type=int, default=1, help="Numero de iteraciones")

  args = parser.parse_args()
  filename = args.file
  iterations = args.iter

  mesh = om.read_polymesh(filename)
  new_mesh = catmull_clark_iter(mesh, iterations)

  ps.init()
  _ = ps.register_surface_mesh(
    "mesh_original", 
    mesh.points(), 
    mesh.face_vertex_indices(),
    transparency=0.5,
    enabled=False
    )
  
  _ = ps.register_point_cloud(
    "point_cloud", 
    new_mesh.points(),
    radius=0.01,
    enabled=False
    )
  
  _ = ps.register_surface_mesh(
    "new_mesh",
    new_mesh.points(),
    new_mesh.face_vertex_indices(),
    edge_width=1,
  )

  ps.show()

if __name__ == "__main__": main()