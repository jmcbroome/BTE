from setuptools import Extension, setup
from Cython.Build import cythonize
import subprocess

#build the tbb library. It comes with a makefile prepared.
subprocess.run("make -j",cwd="./dependencies/tbb/build",shell=True,check=True)
#build the matree protoc reader files for this system and put them in the src
subprocess.run("protoc --cpp_out=./src parsimony.proto",shell=True,check=True)

extensions = [Extension("mat",["src/mat.pyx"],
    include_dirs=["${CONDA_PREFIX}/include/google/protobuf", "${CONDA_PREFIX}/include/", "./dependencies/oneTBB/include"],
    library_dirs=['${CONDA_PREFIX}/lib/','./dependencies/oneTBB/build/'],
    libraries=['tbb', 'protobuf','boost_system','boost_iostreams'],
    language='c++'
    )]

setup(name='mat',ext_modules=cythonize(extensions,language_level=3))