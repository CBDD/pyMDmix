.. pyMDMix documentation master file, created by
   sphinx-quickstart on Thu Feb  6 09:14:42 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Replicas Module documentation
===================================

.. automodule:: Replicas

Auxiliary functions
-------------------
.. autofunction:: Replicas.loadReplica
.. autofunction:: Replicas.createReplicas

Replica class
-------------
.. .. autoclass:: Replicas.Replica(solvent, name, fromfile, off, crd, top, ...) 
.. autoclass:: Replicas.Replica
	:member: __init__()
	:special-members: __init__

Replica attributes
------------------

.. .. autoattribute:: Replicas.Replica.name
.. .. autoattribute:: Replicas.Replica.solvent
.. .. autoattribute:: Replicas.Replica.off 
.. .. autoattribute:: Replicas.Replica.unitname 
.. .. autoattribute:: Replicas.Replica.top 
.. .. autoattribute:: Replicas.Replica.pdb 
.. .. autoattribute:: Replicas.Replica.crd 
.. .. autoattribute:: Replicas.Replica.analyzenetcdf 
.. .. autoattribute:: Replicas.Replica.path 
.. .. autoattribute:: Replicas.Replica.extraResidues 
.. .. autoattribute:: Replicas.Replica.restrMask 
.. .. autoattribute:: Replicas.Replica.alignMask 
.. .. autoattribute:: Replicas.Replica.referencePdb 
.. .. autoattribute:: Replicas.Replica.replFilePath 

Replica methods
---------------

.. automethod:: Replicas.Replica
	:members:
	
Examples
--------

Importing existing data
^^^^^^^^^^^^^^^^^^^^^^^

In this example we will create an empty replica folder and add existing 
sources from an imaginary previous simulation. It was run with ethanol mixture (named *ETA*) for 40ns.

.. code-block:: python

	>>> previousdata = {'pdb':'/some/path/system.pdb', 'crd':'/some/path/system.crd','top':'/some/path/system.top','mdfolder':'/some/path/production','eqfolder':'/some/path/equilibration'}
	>>> from Replicas import Replica
	>>> replica = Replica(name='mynewreplica', nanos=40, solvent='ETA')
	>>> replica.createFolder()
	>>> replica.importData(**previousdata) # This will link all existing files inside the created folders	

Guide to developers
-------------------
Replica flexible configuration (developers guide)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`Replica` objects can be configured by mean of their constructor method.
Arguments not present in at construction time, will take default values from default settings
and user settings **replica-settings.cfg** files.

.. code-block:: python

    >>>from Replicas import Replica
    >>>replica = Replica('customreplica', nanos=40, temp=298)
    >>>print replica.nanos
    40
    >>>print replica.temp
    298
    >>>defreplica = Replica()
    >>>print defreplica.nanos
    20
    >>>print defreplica.temp
    300
    
*temp* and *nanos* attributes where assigned from user configuration file which should be at user's home
directory $HOME/.mdmix/replica-settings.cfg. Values not explicitely assigned there, will be taken from
default configurations in package installation directory.

This system for building instances permits the developer to modify/add/remove attributes to the instance without
modifying any code. For instance, if a new pair *int-myattr=200* is written in user's replica configuration file,
default replicas will also have that new attribute.

