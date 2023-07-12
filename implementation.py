
import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import eigs, spsolve
import openmesh

def calcular_baricentro(mesh, face):
  points = mesh.points()
  mean_vertex = np.array([0, 0, 0], dtype=float)
  vertex_count = 0.0

  #iterador de vertices de la cara
  for v in mesh.fv(face):
    mean_vertex += points[v.idx()]
    vertex_count += 1.0

  mean_vertex /= vertex_count

  return mean_vertex

def calcular_aricentro(mesh, edge, new_mesh):
  points = mesh.points()
  new_mesh_points = new_mesh.points()
  
  he = mesh.halfedge_handle(edge, 0)
  mean_vertex = np.array([0, 0, 0], dtype=float)
  vertex_count = 0.0

  #detectar si la arista no es frontera
  if not mesh.is_boundary(edge):
    #obtener caras adyacentes
    f1 = mesh.face_handle(he)
    f2 = mesh.face_handle(mesh.opposite_halfedge_handle(he))

    #sumar los baricentros de las caras
    mean_vertex += new_mesh_points[f1.idx()] + new_mesh_points[f2.idx()]
    vertex_count += 2.0

  #obtener vertices de la arista
  v1 = mesh.from_vertex_handle(he)
  v2 = mesh.to_vertex_handle(he)
  
  #sumar los vertices de la arista
  mean_vertex += points[v1.idx()] + points[v2.idx()]
  vertex_count += 2.0

  #promediar
  mean_vertex /= vertex_count

  return mean_vertex

def calcular_esquinas(mesh, vertex, new_mesh):
  points = mesh.points()
  new_mesh_points = new_mesh.points()
  
  v_new = np.array([0, 0, 0], dtype=float)
  v_old = points[vertex.idx()]
  edges_count = 0.0
  faces_count = 0.0

  if not mesh.is_boundary(vertex):
    edge_vertex = np.array([0, 0, 0], dtype=float)
    #circulador de aristas
    for eh in mesh.ve(vertex):
      index = eh.idx()
      aricentro = new_mesh_points[mesh.n_faces() + index]
      edge_vertex += (aricentro - v_old)
      edges_count += 1.0

    face_vertex = np.array([0, 0, 0], dtype=float)
    #iterador de caras
    for f in mesh.vf(vertex):
      index = f.idx()
      baricentro = new_mesh_points[index]
      face_vertex += (baricentro - v_old)
      faces_count += 1.0

    edge_vertex /= np.power(edges_count, 2)
    face_vertex /= np.power(edges_count, 2)

    v_new = v_old + edge_vertex + face_vertex

  else:
    edge_vertex = np.array([0, 0, 0], dtype=float)
    #circulador de aristas frontera
    for eh in mesh.ve(vertex):
      if mesh.is_boundary(eh):
        index = eh.idx()
        aricentro = new_mesh_points[mesh.n_faces() + index]
        edge_vertex += (aricentro - v_old)
        edges_count += 1.0

    edge_vertex /= np.power(edges_count, 2)
    v_new = v_old + edge_vertex

  return v_new

def catmull_clark(mesh: openmesh.PolyMesh):
  #En cada iteración:
    #Construir vertices de caras: baricentros
    #Construir vertices de aristas: promedio de los vertices de las aristas
    #Actualizar posiciones de los vertices originales

  #crear una nueva malla
  #indices:
    #[0, ..., n_faces() - 1]: baricentros
    #[n_faces(), ..., n_faces() + n_edges() - 1]: aricentros
    #[n_faces() + n_edges(), ..., n_faces() + n_edges() + n_vertices() - 1]: vertices originales
  new_mesh = openmesh.PolyMesh()
  n_faces = mesh.n_faces()
  n_edges = mesh.n_edges()

  #iterador de caras
  for f in mesh.faces():
    #agregar el vertice a la nueva malla
    _ = new_mesh.add_vertex(calcular_baricentro(mesh, f))

  #iterador de aristas
  for e in mesh.edges():
    #agregar el vertice a la nueva malla
    _ = new_mesh.add_vertex(calcular_aricentro(mesh, e, new_mesh))

  #iterador de vertices
  for v in mesh.vertices():
    #agregar el vertice a la nueva malla
    _ = new_mesh.add_vertex(calcular_esquinas(mesh, v, new_mesh))

  for f in mesh.faces():
    baricentro = new_mesh.vertex_handle(f.idx())
    
    vertices = [vh for vh in mesh.fv(f)]
    edges = [eh for eh in mesh.fe(f)]

    # [WARNING] this is supposed to work mostly for quads, some weird stuff happens with triangles and borders sometimes
    # EXAMPLE model airplane_0627 gets some weird holes which messes up the mesh for the following iterations
    # but we'll let it slide for now using the range with len(edges) instead of 4
    for i in range(len(edges)):
      esquina = new_mesh.vertex_handle(
        n_faces + n_edges + vertices[i].idx()
      )
      aricentro_anterior = new_mesh.vertex_handle(
        n_faces + edges[i].idx()
      )
      aricentro_siguiente = new_mesh.vertex_handle(
        n_faces + edges[(i + 1) % len(edges)].idx()
      )

      _ = new_mesh.add_face([esquina, aricentro_anterior, baricentro, aricentro_siguiente])

  return new_mesh

def catmull_clark_iter(mesh: openmesh.PolyMesh, iterations: int):
  new_mesh = mesh
  for _ in range(iterations):
    new_mesh = catmull_clark(new_mesh)

  return new_mesh