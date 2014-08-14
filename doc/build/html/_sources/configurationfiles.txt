.. _pcf:

################################
Project Configuration File (PCF)
################################
File used to create Project instances. Here, we will give information about the system
we are going to simulate. It has this shape:

.. literalinclude:: project.cfg
	:language: ini 

.. _rcf:

################################
Replica Configuration File (RCF)
################################
With same idea, replicas can also be created using a configuration file.
In this case, the majority of information we will give is regarding
the simulation parameters and conditions.

.. literalinclude:: replica.cfg
	:language: ini

.. _completecf:

######################
Complete configuration 
######################
In this case, the file contains both sections. All replicas will take system information from the system and therefore, do not require any system information.

.. literalinclude:: complete.cfg
	:language: ini
