
from distutils.core import setup, Extension
import numpy as np

from Cython.Build import cythonize


extensions = [
    Extension("ellipse", ["ellipse.pyx"]),
    Extension("time_connected_clusters", ["time_connected_clusters.pyx"]),
]

setup(
    include_dirs = [np.get_include()],
    ext_modules = cythonize(extensions, include_path=[np.get_include()]),
)
