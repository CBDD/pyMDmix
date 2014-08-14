.. pyMDMix documentation master file, created by
   sphinx-quickstart on Thu Feb  6 09:14:42 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _settings:

Setting General Options and Attributes through Configuration files
==================================================================

Several classes automatically take arguments from default configuration files distributed along the package.
Once pyMDMix is started for the first time, it will also make a copy of these configuration files inside the user's home directory
for easy modification. If any parameter is modified by the user, it will have higher priority and the default one will be ignored.
For restoring initial file, just remove it from the user directory.

Two configuration files govern the program::
	- General settings (``settings.cfg``)
	- MD settings (``md-settings.cfg``)

General settings
----------------

This is the default file for configuring general and project options in pyMDMix.
It can be found at the package installation directory (``$INSTALLDIR``) under ``$INSTALLDIR/data/defaults/settings.cfg``
or at user's home directory ``~./mdmix/settings.cfg``

.. literalinclude:: settings.cfg
	:language: ini 

.. _mdsettings:

MD default settings
------------------------
This settings file contains two kind of information::
	- MD parameters: control default integration timestep, number of nanoseconds, temperature, trajectory writing frequency or default restraining schema
	- Folder and replica control parameters: names for the folders or output filename templates.

The current built-in settings file contains::
	.. literalinclude:: md-settings.cfg
		:language: ini

.. _settingsmodule:

Settings Module
===============

.. automodule:: settings
