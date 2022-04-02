# BTE
Cython API for the [Mutation Annotated Tree (MAT) Online Phylogenetics Toolkit](https://github.com/yatisht/usher) (In Early Development). Developer version distributed as source. 

# Overview
This repository will allow the user to leverage the power of the Mutation Annotated Tree file format and library in their Python scripts, allowing for efficient and effective analysis of global SARS-CoV-2 and other pathogen phylogenies. 

UCSC maintains a daily-updated mutation-annotated tree in protobuf format containing all publicly available SARS-CoV-2 genome sequences [here](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/). To try out this tool, download the latest tree, build the library, and jump straight to your Python analysis!

## Quickstart
If you're on macOS and have conda installed, you can install our package from conda with the following:

```
conda install -c jmcbroome bte
```

And proceed directly to your analysis by calling "import mat" in your python scripts!

We intend to add more platforms and python versions soon!

## Build From Source Instructions

You may need to build this library from source if you are adding functionality or on a currently unsupported architecture.

### First Time Setup

This project is dependent on a few key libraries that need to be available for linking. We use conda for environment management, so you will need miniconda or anaconda. 

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
```

If the .yml isn't working for you, you can try:

```
conda create --name bte
conda activate bte
conda install -c conda-forge -c anaconda protobuf boost-cpp cython tbb-devel=2019.0
```

When building your own environment, you will also need the general suite of C++ compiler tools. On a mac, you may need xcode CLI. On linux, you may need other compiler tools. You may need to call 

```
conda install -c conda-forge cxx-compiler
conda install -c anaconda make 
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

From a script in the same directory as the .so!

### Example Usage

We provide a unit test script you can use to validate that the library is functioning as intended.

```
python3 -m unittest run_test.py
```

Additionally, you can explore a simple example analysis in the BTE Tutorial jupyter notebook included in this repository! Just
retrieve the latest tree from the UCSC repository before you being.

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
```