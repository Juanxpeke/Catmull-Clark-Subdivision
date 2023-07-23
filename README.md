# Catmull-Clark-Subdivision
Project for geometric processing course. It shows a simple implementation of the Catmull-Clark subdivision algorithm using OpenmMesh for Python.

## Build
You just need to clone this repository and a Python environment with openmesh and polyscope modules. For this, Python 3.8 version is recommended, as it works perfectly with both libraries.

## Run

For a default run just type in the terminal.

```python visualization.py```

Otherwise, you can specify the model and the level of detail given by the number of subdivisions applied to the model.

```python visualization.py --file filename --iter number_of_iterations```
