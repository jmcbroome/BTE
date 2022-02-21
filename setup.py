from setuptools import Extension, setup
from Cython.Build import cythonize
import subprocess
from os import environ

CONDA_PREFIX = environ.get('CONDA_PREFIX', None)
if CONDA_PREFIX == None:
    print("Can't find conda installation- aborting.")
    exit(1)

#build the matree protoc reader files for this system and put them in the src
subprocess.run("protoc --cpp_out=./src parsimony.proto",shell=True,check=True)

extensions = [Extension("mat",["src/mat.pyx"],
    include_dirs=[CONDA_PREFIX+"/include/google/protobuf", CONDA_PREFIX+"/include/"],
    library_dirs=[CONDA_PREFIX+'/lib/'],
    libraries=['tbb', 'protobuf','boost_system','boost_iostreams'],
    language='c++'
    )]  

setup(name='mat',ext_modules=cythonize(extensions,language_level=3))