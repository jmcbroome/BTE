FROM continuumio/miniconda3:22.11.1
USER root
WORKDIR /home
RUN git clone --recurse-submodules https://github.com/jmcbroome/BTE.git
WORKDIR BTE
RUN conda env create -f bte.yml
SHELL ["conda", "run", "-n", "bte", "/bin/bash", "-c"]
RUN python3 setup.py build_ext && \
 python3 setup.py install && \
 ln -s /home/BTE/build/lib.linux-x86_64-cpython-310/bte.cpython-310-x86_64-linux-gnu.so /opt/conda/lib/python3.10/site-packages