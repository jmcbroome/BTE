# matreePy
Cython API for the Mutation Annotated Tree (MAT) Online Phylogenetics Toolkit (In Early Development)

## Dependencies
Boost

Protobuf

tbb 2019_U6 (included here)


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

In the cloned directory before proceeding. Create a conda environment containing our other key dependencies.
```
conda env create -f matreepy.yml
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

