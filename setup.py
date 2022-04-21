from setuptools import Extension, setup
from Cython.Build import cythonize
import subprocess
from os import environ

CONDA_PREFIX = environ.get('CONDA_PREFIX', None)
if CONDA_PREFIX == None:
    print("Can't find conda installation- aborting.")
    exit(1)

#build the protoc reader files for this system and put them in the src
subprocess.run("protoc --cpp_out=../ parsimony.proto",shell=True,check=True,cwd='src/usher')

extensions = [Extension("bte",["src/bte.pyx"],
    include_dirs=[CONDA_PREFIX+"/include/google/protobuf", CONDA_PREFIX+"/include/"],
    library_dirs=[CONDA_PREFIX+'/lib/'],
    libraries=['tbb', 'protobuf','boost_system','boost_iostreams'],
    language='c++',
    extra_compile_args=["-std=c++17"]
    )]

setup(name='bte',ext_modules=cythonize(extensions,language_level=3,annotate=True))
