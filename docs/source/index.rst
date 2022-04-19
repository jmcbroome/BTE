.. BTE documentation master file, created by
   sphinx-quickstart on Thu Apr 14 17:09:13 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the BTE documentation!
=================================

BTE (Big Tree Explorer) is a Python extension for analysis and traversal of extremely large phylogenetic trees. It's based on the
Mutation Annotated Tree (MAT) library which underlies UShER and matUtils, premier software for pandemic-scale phylogenetics. 

This tool is generally intended as a replacement for ETE3, Biopython.Phylo, and similar Python phylogenetics packages for Mutation Annotated Trees (MATs). 
Using standard packages with MATs requires conversion to newick and the maintenance of mutation annotations as a separate data structure, 
generally causing inconvenience and slowing both development and runtime. BTE streamlines this process by exposing the heavily optimized MAT library underlying UShER and matUtils to Python, allowing for efficient and convenient use of MATs in a Python development environment!

UCSC maintains a repository, updated each day, containing the complete and latest publicly-available global SARS-CoV-2 phylogenetic tree in MAT protobuf format here. 

http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/

To try out this tool, download the latest tree, build or install the extension, and jump straight to your Python analysis!

Installation
============

On Mac or Linux, you can install BTE by running the following command:

      conda install -c jmcbroome bte

See the README if you encounter installation difficulties or need to build a local version of the extension.

Basic Logic
===========

BTE contains two primary classes, the MATree and the MATNode. 

The MATree class represents a single mutation-annotated phylogenetic tree. It can be created from a MAT protobuf (.pb) file, from a newick and a vcf together, or from an Auspice-format JSON like that used by Nextstrain.
The MATree class has many useful functions for manipulation, subsetting, traversal, and summarization of the tree it represents.

The MATNode class represents a single node member of the Mutation Annotated Tree. It has parent, child, mutation, and annotation attributes
which can be used to traverse or compute statistics from the tree.

The full API documentation follows.

BTE API Documentation
=====================

.. toctree::

:ref:`genindex`

.. automodule:: bte
   :members: