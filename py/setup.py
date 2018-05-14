
from distutils.core import setup, Extension
import numpy as np

from Cython.Build import cythonize


extensions = [
    Extension("ellipse", ["ellipse.pyx"]),
]

setup(
    include_dirs = [np.get_include()],
    ext_modules = cythonize(extensions, include_path=[np.get_include()]),
)
