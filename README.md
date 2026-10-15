## pyMDMix --- http://mdmix.sourceforge.net

The program is distributed under GNU GPLv3 license. Find the license file
under Licenses/ folder.

DOCUMENTATION
=============
All documentation on program usage is online at
http://mdmix.sourceforge.net

INSTALLATION
============

1- Dependencies
---------------
This version of pyMDMix depends on:
  - ambertools>=12
  - python>=3.9
  - biskit>=3.0 (installed automatically; this is the Python 3 release of Biskit, imported as `biskit`)

make sure ambertools environment variables are set

2 - Installation process
------------------------
it is advised to install pyMDMix in a virtual environment

there are three recommended ways to install pyMDMix:
1. from the repository by
`python -m pip install [insert address here]`

2. from the project's local folder after cloning the repo by
`python -m pip install .`

3. Use conda or mamba: we will first create the correct conda environment which will already contain ambertools. We will then set the AMBERHOME variable within the environment and finally proceed to install pymdmix from the local git cloned files. Make sure `$CONDA_PREFIX` points to the conda installation directory (it may happen when you activate the new environment this varaible is lost, in that case change the variable for the explicit path). 

```bash
git clone https://github.com/CBDD/pyMDmix.git
cd pyMDmix
conda env create -f environment.yml
conda activate pymdmix3
conda env config vars set AMBERHOME=$CONDA_PREFIX
conda deactivate && conda activate pymdmix3
pip install -e .
```

3 - Testing
-----------

You can test the package has been successfuly installed by activating
the environment and, from a folder other than the cloned repository, running
`python -m pyMDMix`
if everything went fine, you will see something like
`Welcome to MDMix 2.6`

The main script should have been installed correclty also. check it by running
`mdmix -h`

Test the program works correctly:
	Move again to the package source directory and type:
		> python pyMDMix/test.py all
	This will run a series of source code checks.
	Some tests need a working AMBER installation (tleap, cpptraj, ambpdb) and
	internet access (PDB2PQR web server), so they will fail without them.

Tutorial
--------
A short Spanish tutorial (cajas de solvente OPC, creacion de proyectos, analisis) is in TUTORIAL.md. The step-by-step guide to build mixed-solvent boxes with OPC water (build_box.sh, boxcalc.py) is in TUTORIAL_CAJAS.md.
A complete project configuration template is in pyMDMix/data/templates/complete.cfg
(`mdmix create template -f proyecto.cfg` copies it to the current folder).

Python 3 notes
--------------
pyMDMix runs on Python 3 (>=3.9). Projects, replicas, systems, hot spot sets and the
solvent database saved with the Python 2 version can still be opened: they are read with
the `latin1` encoding and saved again in the Python 3 pickle format the next time they
are written. Once saved by Python 3 they can no longer be read by the Python 2 version.

4 - Using Docker
----------------
Build and run from a docker container using the provided Dockerfile.

- docker build -t pymdmix .
- docker run -v $PWD:/mnt pymdmix -h

The docker run command will mount the current working directory (windows users should replace $PWD by %cd%) so the container has access to the current location. Output will be also in the current working directory. 
s

4 - Enjoy!
----------
Read program usage at online documentation.


