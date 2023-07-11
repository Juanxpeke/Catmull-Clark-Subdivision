
import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import eigs, spsolve
import openmesh
import argparse
import os
import polyscope as ps

def calcular_baricentro(mesh, face):
  points = mesh.points()
  mean_vertex = np.array([0, 0, 0], dtype=float)
  vertex_count = 0
  #iterador de vertices de la cara
  for v in mesh.fv(face):
    mean_vertex += points[v.idx()]
    vertex_count += 1

  mean_vertex /= vertex_count

  return mean_vertex

def calcular_aricentro(mesh, edge, new_mesh):
  points = mesh.points()
  new_mesh_points = new_mesh.points()
  mean_vertex = np.array([0, 0, 0], dtype=float)
  vertex_count = 0.0

  he = mesh.halfedge_handle(edge, 0)
  # detectar si la arista es no es frontera
  if not mesh.is_boundary(edge):
    #obtener caras adyacentes
    f1 = mesh.face_handle(he)
    f2 = mesh.face_handle(mesh.opposite_halfedge_handle(he))

    # sumar los baricentros de las caras
    mean_vertex += new_mesh_points[f1.idx()] + new_mesh_points[f2.idx()]
    vertex_count += 2.0

  #obtener vertices de la arista
  v1 = mesh.from_vertex_handle(he)
  v2 = mesh.to_vertex_handle(he)
  
  #sumar los vertices de la arista
  mean_vertex += points[v1.idx()] + points[v2.idx()]
  vertex_count += 2.0

  # promediar
  mean_vertex /= vertex_count

  return mean_vertex

def calcular_esquinas(mesh, vertex, new_mesh):
  points = mesh.points()
  new_mesh_points = new_mesh.points()
  
  v_new = np.array([0, 0, 0], dtype=float)
  v_old = points[vertex.idx()]
  edges_count = 0
  faces_count = 0

  edge_vertex = np.array([0, 0, 0], dtype=float)
  #circulador de aristas
  for eh in mesh.ve(vertex):
    index = eh.idx()
    aricentro = new_mesh_points[mesh.n_faces() + index]
    edge_vertex += (aricentro - v_old)
    edges_count += 1
  edge_vertex /= np.power(edges_count, 2)

  face_vertex = np.array([0, 0, 0], dtype=float)
  #iterador de caras
  for f in mesh.vf(vertex):
    index = f.idx()
    baricentro = new_mesh_points[index]
    face_vertex += (baricentro - v_old)
    faces_count += 1
  
  if not mesh.is_boundary(vertex):
    face_vertex /= np.power(edges_count, 2)
    v_new = (v_old + edge_vertex + face_vertex)
  else:
    v_new = v_old

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
    esquina_handle = new_mesh.add_vertex(calcular_esquinas(mesh, v, new_mesh))

    #aristas
    first_handle = None
    edge_count = 0
    for eh in mesh.ve(v):
      edge_count += 1
      if first_handle is None:
        first_handle = eh
        last_handle = eh
        continue
      
      if last_handle != eh and mesh.is_boundary(last_handle):
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

      baricentro_common_handle = None

      if f1.idx() == -1:
        baricentro_common_handle = new_mesh.vertex_handle(f2.idx())
      elif f2.idx() == -1:
        baricentro_common_handle = new_mesh.vertex_handle(f1.idx())
      else:
        baricentro_f1_handle = new_mesh.vertex_handle(f1.idx())
        baricentro_f2_handle = new_mesh.vertex_handle(f2.idx())

        a1 = new_mesh.vertex_handle(last_f1.idx())
        a2 = new_mesh.vertex_handle(last_f2.idx())
        if a1 == baricentro_f1_handle or a1 == baricentro_f2_handle:
          baricentro_common_handle = a1
        else:
          baricentro_common_handle = a2

      #agregar cara [WIP]
      
      _ = new_mesh.add_face([esquina_handle, aricentro_handle, baricentro_common_handle, last_aricentro_handle])
    
      # [WIP]
      last_handle = eh

    if edge_count == 0:
      continue

    if mesh.is_boundary(last_handle):
      continue

    # falta agregar la ultima cara con la primera arista
    eh_index = first_handle.idx()
    last_index = last_handle.idx()
    he_handle = mesh.halfedge_handle(first_handle, 0)
    last_he_handle = mesh.halfedge_handle(last_handle, 0)

    aricentro_handle = new_mesh.vertex_handle(mesh.n_faces() + eh_index)
    last_aricentro_handle = new_mesh.vertex_handle(mesh.n_faces() + last_index)

    #caras adyacentes a la arista
    f1 = mesh.face_handle(he_handle)
    f2 = mesh.face_handle(mesh.opposite_halfedge_handle(he_handle))
    #caras adyacentes a la arista anterior
    last_f1 = mesh.face_handle(last_he_handle)
    last_f2 = mesh.face_handle(mesh.opposite_halfedge_handle(last_he_handle))

    baricentro_common_handle = None

    if f1.idx() == -1:
      baricentro_common_handle = new_mesh.vertex_handle(f2.idx())
    elif f2.idx() == -1:
      baricentro_common_handle = new_mesh.vertex_handle(f1.idx())
    else:
      baricentro_f1_handle = new_mesh.vertex_handle(f1.idx())
      baricentro_f2_handle = new_mesh.vertex_handle(f2.idx())

      a1 = new_mesh.vertex_handle(last_f1.idx())
      a2 = new_mesh.vertex_handle(last_f2.idx())
      if a1 == baricentro_f1_handle or a1 == baricentro_f2_handle:
        baricentro_common_handle = a1
      else:
        baricentro_common_handle = a2

    #agregar cara [WIP]
    _ = new_mesh.add_face([esquina_handle, aricentro_handle, baricentro_common_handle, last_aricentro_handle])

  return new_mesh

def catmull_clark_iter(mesh: openmesh.PolyMesh, iterations: int):
  new_mesh = mesh
  for i in range(iterations):
    new_mesh = catmull_clark(new_mesh)

  return new_mesh
  

def catmull_clark2(mesh: openmesh.PolyMesh):
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
    #iterador de vertices de la cara
    for v in mesh.fv(f):
      esquina = new_mesh.vertex_handle(
        n_faces + n_edges + v.idx()
      )
      
      # [TODO]
      #obtener arista anterior y siguiente
      #obtener aricentro anterior y siguiente
      
      #agregar cara
      # [TODO]

  return new_mesh

def catmull_clark_iter2(mesh: openmesh.PolyMesh, iterations: int):
  new_mesh = mesh
  for i in range(iterations):
    new_mesh = catmull_clark(new_mesh)

  return new_mesh