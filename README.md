[![Documentation Status](https://readthedocs.org/projects/shesha/badge/?version=py3)](http://shesha.readthedocs.io/en/py3/?badge=py3) [![Anaconda-Server Badge](https://anaconda.org/compass/compass/badges/installer/conda.svg)](https://conda.anaconda.org/compass)

Table of Contents
=================

  * [Requirements](#requirements)
  * [Installation de COMPASS via conda](#installation-de-compass-via-conda)
  * [Installation de SHESHA package for COMPASS](#installation-de-shesha-package-for-compass)
  * [More documentation (maybe not fully up-to-date)](#more-documentation-maybe-not-fully-up-to-date)
  * [Questions?](#questions)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)

## Requirements

Linux computer with CUDA 9.0

## Installation de COMPASS via conda

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
export PATH=/home/sevin/miniconda3/bin:$PATH
conda install -c compass compass -y
```

## Installation de SHESHA package for COMPASS

```bash
cd
git clone https://github.com/ANR-COMPASS/shesha.git
export SHESHA_ROOT=$HOME/shesha
export PYTHONPATH=$SHESHA_ROOT/src:$PYTHONPATH
export PYTHONDONTWRITEBYTECODE=1
cd $SHESHA_ROOT
ipython -i test/closed_loop.py data/par/par4bench/scao_sh_16x16_8pix.py
```

## More documentation (maybe not fully up-to-date)

doc auto-generated from code: http://shesha.readthedocs.io

wiki page of the COMPASS project: https://projets-lesia.obspm.fr/projects/compass/wiki/Wiki

## Questions?

Please feel free to create an issue on Github for any questions and inquiries.
