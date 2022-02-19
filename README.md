# matreePy
Cython API for the Mutation Annotated Tree (MAT) Online Phylogenetics Toolkit (In Early Development)

## Dependencies
Boost

Protobuf

tbb 2019_U6 (included here)


## Build Instructions

### First Time Setup

This project is dependent on a few key libraries that need to be available for linking. oneTBB is included as a submodule given our dependency on an older version of the library. To ensure oneTBB is available to be build, when cloned, use

```
git clone --recurse-submodules https://github.com/jmcbroome/matreePy
```

If you forget, you can use 

```
git submodule update --init
```

Before proceeding to initialize the oneTBB source submodule.

```
python3 setup.py build_ext --inplace
python3 test.py
```
