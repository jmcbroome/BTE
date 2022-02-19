# matreePy
Cython API for the Mutation Annotated Tree (MAT) Online Phylogenetics Toolkit (In Early Development). Developer version distributed as source. 

## Core Dependencies

Cython 

Boost

Protobuf

tbb 2019_U6 (included here)

Conda, while not strictly required, is very highly recommended for environment management.

You will also need the general suite of C++ compiler tools. On a mac, you will need xcode CLI. On linux, you may need other compiler tools. You may want to call 
```
conda install -c conda-forge cxx-compiler
conda install -c anaconda make 
```
If you don't have make/cmake, clang/gcc and related available.

## Build Instructions

### First Time Setup

This project is dependent on a few key libraries that need to be available for linking. We use conda for environment management, so you will need miniconda or anaconda. oneTBB is included as a submodule given our dependency on an older version of the library. To ensure oneTBB is available to be built, when cloned, use

```
git clone --recurse-submodules https://github.com/jmcbroome/matreePy
```

If you forget, you can use 

```
git submodule update --init
```

After the submodule is available, it needs to be built for your system. Navigate to its directory and apply make.

```
cd dependencies/oneTBB/
make -j
cd ../..
```

If this step fails because of "missing architecture" when building tbb, you can try manually setting the architecture to use- add "arch=x86_64" to the make call for example.

The rest of our key dependencies are handled by conda.

```
conda env create -f matreepy.yml
```

If the .yml isn't working for you (which it may not at this time), you can try:

```
conda create --name matreepy
conda activate matreepy
conda install -c anaconda protobuf
conda install -c conda-forge boost-cpp cython
```

### Building the Python-importable library

Once all libraries are available, proceed to compile the .so.
```
python3 setup.py build_ext --inplace
```

If successful, you can now import from the .so into your python environment with 

```
import mat
```

From a script in the same directory as the .so! Some basic syntax examples are available in the test.py; further updates to the API for usability and tutorials are forthcoming.

