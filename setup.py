from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [Extension("mat",["mat.pyx"],
    include_dirs=["/Users/jmcbr/opt/anaconda3/include/google/protobuf", "/opt/homebrew/include", "./build/oneTBB/include",'/opt/homebrew/Cellar/boost/1.76.0/include/'],
    library_dirs=['./build/oneTBB/build/','/opt/homebrew/Cellar/boost/1.76.0/lib/','/opt/homebrew/opt/boost/lib/','/Users/jmcbr/opt/anaconda3/envs/usher/lib/'],
    libraries=['tbb', 'protobuf','boost_system','boost_iostreams'],
    #extra_compile_args=['-std=c++11'],
    language='c++'
    )]

setup(name='mat',ext_modules=cythonize(extensions,language_level=3))