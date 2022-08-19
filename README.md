# BIG TREE EXPLORER (BTE)
<img src="bte_logo.png" width="500" height="500">

Cython API for the [Mutation Annotated Tree (MAT) Online Phylogenetics Toolkit](https://github.com/yatisht/usher). 

# Overview
This package will allow the user to leverage the power of the Mutation Annotated Tree file format and library in their Python scripts, allowing for efficient and effective analysis of global SARS-CoV-2 and other pathogen phylogenies. 

This package is generally intended as a replacement for ETE3, Biopython.Phylo, and similar Python phylogenetics packages for Mutation Annotated Trees (MATs). Using standard packages with MATs requires conversion to newick and the maintenance of mutation annotations as a separate data structure, generally causing inconvenience and slowing both development and runtime. BTE streamlines this process by exposing the heavily optimized MAT library underlying UShER and matUtils to Python, allowing for efficient and convenient use of MATs in a Python development environment!

UCSC maintains a repository, updated each day, containing the complete and latest publicly-available global SARS-CoV-2 phylogenetic tree in MAT protobuf format [here](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/). To try out this tool, download the latest tree, build the library, and jump straight to your Python analysis!

We also provide a [binder](https://github.com/jmcbroome/bte-binder) of example analyses, which you can launch by clicking on the badge below!

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jmcbroome/bte-binder/HEAD)

We also provide complete [documentation](https://jmcbroome.github.io/BTE/build/html/index.html) for all available functions!

## Who is this for?

If you are 

- a researcher who wants to query the global SARS-CoV-2 phylogeny for your scientific work in a way that is more complex than what is directly supported by UShER and matUtils
- an epidemiologist or public health officer who wants to perform a hands-on Python analysis examining genetic data specific to your area or time of interest
- a developer who wants the convenience and fast development time of Python without having to deal with the overhead and frustration of file conversions and newick manipulation

This package is for you!

## Quickstart

If you're on osx64 (MacOS) or linux64 (most Linux distributions) and have conda installed, you can install our package via the bioconda channel.

```
conda install -c bioconda bte
```

Download the latest public SARS-CoV-2 tree:

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
```

And proceed directly to your analysis in Python!

```
import bte
tree = bte.MATree("public-latest.all.masked.pb.gz")
```

### A Note on Versions and Architectures

We provide conda builds for Linux and MacOS with Python versions >=3.8. If you're on Windows 10+, you can install BTE on the [Linux subsystem](https://docs.microsoft.com/en-us/windows/wsl/about) and perform analyses and run notebooks from the subsystem. If you're on an earlier version of Python and unable to update, you may still be able to build a local extension using the instructions below.

## Build From Source Instructions

You may need to build this package from source if you are adding functionality, on a currently unsupported architecture, on an unsupported version of Python, or just if you want the latest and greatest version. 

### Set Up Build Environment

This package is dependent on a few key libraries that need to be available for linking. We use conda for environment management, so you will need miniconda or anaconda. 

The first step is to clone this repository. This repository relies on UShER as a submodule dependency, so --recurse-submodules
needs to be added to the call.

```
git clone --recurse-submodules https://github.com/jmcbroome/BTE
```

If you forgot to add this argument to your git clone call, you can run the following in the cloned repository.

```
git submodule update --init
```

The next step is to prepare the environment and all relevant dependencies.

```
conda env create -f bte.yml
conda activate bte
```

### Compile and Build Module

Once all libraries are available, proceed to compile the shared object file (.so).

```
python3 setup.py build_ext --inplace
```

If successful, you can now import from the .so into your python environment with 

```
import bte
```

From a script in the same directory as the .so!

We provide a unit test script you can use to validate that the package is functioning as intended.

```
python3 -m unittest run_test.py
```

## Installation Issues

### Installation from `conda` Channel

BTE, as a python extension, is version specific. If you're unable to "import bte" after successful installation via conda, ensure you're using the version of Python associated with your conda installer (python --version). If you encounter an error related to missing ".so" files, setup.py build_ext raises exceptions with Protoc versions, or conda is unable to resolve environment conflicts, you may have a broken environment. The fastest solution to any of these related issues is to build a fresh environment with a supported version of Python and install BTE there.

```
conda create --name bte -c conda-forge -c bioconda bte python=3.8
```

### Installation by Local Build

If you're having trouble compiling a local extension (using setup.py), you need the general suite of C++ compiler tools. On a mac, you may need xcode CLI. On linux, you may need other compiler tools. You may need to call 

```
conda install -c conda-forge -c anaconda cxx-compiler make
```
