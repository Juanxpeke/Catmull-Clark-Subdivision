
import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import eigs, spsolve
import openmesh
import argparse
import os
import polyscope as ps

def calcular_baricentro(mesh, face):
  points = mesh.points().copy()
  mean_vertex = None
  vertex_count = 0
  #iterador de vertices de la cara
  for v in mesh.fv(face):
    if mean_vertex is None: 
      mean_vertex = points[v.idx()]
    else: 
      mean_vertex += points[v.idx()]
    vertex_count += 1
    
  if vertex_count < 3: 
    raise Exception(f"Face with less than 3 vertices ({vertex_count})")

  mean_vertex /= vertex_count

  return mean_vertex

def calcular_aricentro(mesh, edge, new_mesh_points):
  points = mesh.points()
  he = mesh.halfedge_handle(edge, 0)

  #obtener caras adyacentes
  f1 = mesh.face_handle(he)
  f2 = mesh.face_handle(mesh.opposite_halfedge_handle(he))

  mean_vertex = new_mesh_points[f1.idx()] + new_mesh_points[f2.idx()]

  #obtener vertices de la arista
  v1 = mesh.from_vertex_handle(he)
  v2 = mesh.to_vertex_handle(he)
  
  mean_vertex = mean_vertex + points[v1.idx()] + points[v2.idx()]
  mean_vertex /= 4

  return mean_vertex

def calcular_esquinas(mesh, vertex, new_mesh_points):
  points = mesh.points()
  v_old = points[vertex.idx()]
  edges_count = 0

  edge_vertex = np.array([0, 0, 0], dtype=float)
  #circulador de aristas
  for eh in mesh.ve(vertex):
    index = eh.idx()
    aricentro = new_mesh_points[mesh.n_faces() + index]
    edge_vertex += (aricentro - v_old)
    edges_count += 1

  edge_vertex /= (edges_count * edges_count)
  
  face_vertex = np.array([0, 0, 0], dtype=float)
  #iterador de caras
  for f in mesh.vf(vertex):
    index = f.idx()
    baricentro = new_mesh_points[index]
    face_vertex += (baricentro - v_old)
  
  face_vertex /= (edges_count * edges_count)

  return v_old + edge_vertex + face_vertex


def catmull_clark(mesh):
  #En cada iteración:
    #Construir vertices de caras: baricentros
    #Construir vertices de aristas: promedio de los vertices de las aristas
    #Actualizar posiciones de los vertices originales

  #crear una nueva malla
  #indices:
    #[0, ..., n_faces() - 1]: baricentros
    #[n_faces(), ..., n_faces() + n_edges() - 1]: aricentros
    #[n_faces() + n_edges(), ..., n_faces() + n_edges() + n_vertices() - 1]: vertices originales
  new_mesh = openmesh.TriMesh()

  #iterador de caras
  for f in mesh.faces():
    mean_vertex = calcular_baricentro(mesh, f)
    #agregar el vertice a la nueva malla
    _ = new_mesh.add_vertex(mean_vertex)

  new_mesh_points = new_mesh.points()
  edge_index = 0
  #iterador de aristas
  for e in mesh.edges():
    mean_vertex = calcular_aricentro(mesh, e, new_mesh_points)
    #agregar el vertice a la nueva malla
    _ = new_mesh.add_vertex(mean_vertex)
    edge_index += 1

  new_mesh_points = new_mesh.points()
  #iterador de vertices
  for v in mesh.vertices():
    esquina_vertex = calcular_esquinas(mesh, v, new_mesh_points)

    #agregar el vertice a la nueva malla
    esquina_handle = new_mesh.add_vertex(esquina_vertex)

    #aristas
    first_handle = None
    for eh in mesh.ve(v):
      if first_handle is None:
        first_handle = eh
        last_handle = eh
        continue

      eh_index = eh.idx()
      last_index = last_handle.idx()
      he_handle = mesh.halfedge_handle(eh, 0)
      last_he_handle = mesh.halfedge_handle(last_handle, 0)

      aricentro_handle = new_mesh.vertex_handle(mesh.n_faces() + eh_index)
      last_aricentro_handle = new_mesh.vertex_handle(mesh.n_faces() + last_index)

      #caras adyacentes a la arista
      f1 = mesh.face_handle(he_handle)
      f2 = mesh.face_handle(mesh.opposite_halfedge_handle(he_handle))
      #caras adyacentes a la arista anterior
      last_f1 = mesh.face_handle(last_he_handle)
      last_f2 = mesh.face_handle(mesh.opposite_halfedge_handle(last_he_handle))

      baricentro_f1_handle = new_mesh.vertex_handle(f1.idx())
      baricentro_f2_handle = new_mesh.vertex_handle(f2.idx())

      baricentro_common_handle = None
      a1 = new_mesh.vertex_handle(last_f1.idx())
      a2 = new_mesh.vertex_handle(last_f2.idx())
      if a1 != baricentro_f1_handle and a1 != baricentro_f2_handle:
        baricentro_common_handle = a2
      else:
        baricentro_common_handle = a1

      #agregar cara [WIP]
      
      _ = new_mesh.add_face([esquina_handle, aricentro_handle, baricentro_common_handle, last_aricentro_handle])
    
      # [WIP]
      last_handle = eh

    # falta agregar la ultima cara con la primera arista

  return new_mesh
