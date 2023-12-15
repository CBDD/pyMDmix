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
  - python>=2.7

make sure ambertools environment variables are set

2 - Installation process
------------------------
it is advised to install pyMDMix in a virtual environment

there are three recommended ways to install pyMDMix:
1. from the repository by
`python -m pip install [insert address here]`

2. from the project's local folder after cloning the repo by
`python -m pip install .`

3. Use conda or mamba: `conda env create -f environment_p27.yml` then activate it `conda activate mdmix-env`

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
	No test should fail.

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


