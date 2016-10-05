Installing the Cubit/Trelis Plugin
==================================

There is a common location that store all the plugins that are available for use. Go 
to `DAGMC Plugins <http://go.wisc.edu/dagmc-trelis>`_ page and download the plugin appropriate
for your operating system.

Linux Install Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Obtain the Linux pluging from the link above, the following operations should be 
performed as the root user (sudo).

1. unpack tar ball in /opt/Trelis-16.0/bin/plugins
2. change directory to /opt/Trelis-16.0/bin/plugins/dagmc
3. Run the install.sh script: ./install.sh

OS/X Install Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

Obtain the OS/X plugin from the link above, the following operations should be 
performed as the root user (sudo).

1. If it does not already exist, create a folder called "plugins" in
   /Applications/Trelis-16.0/Contents/MacOS/
2. unpack tar ball in /Applications/Trelis-16.0/Contents/MacOS/plugins
3. change directory to /Applications/Trelis-16.0/Contents/MacOS/plugins/dagmc
4. Run the install.sh script: ./install.sh

Windows Install Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A plugin for Windows is currently being developed, it will be posted to the 
above link when complete.

Using the Plugin
~~~~~~~~~~~~~~~~

Load your geometry into Cubit/Trelis following the instructions in 
`geometry preparation <workflow/cubit_trelis_workflow.html>`_, markup your
geometry according to the workflow of your choice, for example 
`UWUW  <workflow/uw2.html>`_, and Imprint and Merge. You can now export the 
geometry, ready for use in simulations using
::

   CUBIT> export dagmc <filename> [faceting_tolerance <faceting tolerance>]
                                  [length_tolerance <length tolerance>]
                                  [normal_tolerance <normal tolerance>]
                                  [verbose] [fatal_on_curves]

Where faceting_tolerance, normal_tolerance, and length tolerance are optional arguments. For example, export
the currently loaded file to a file called "large_facets.h5m" with a faceting tolerance of 1.e-5
::

   CUBIT> export dagmc "large_facets.h5m" faceting_tolerance 1.e-5

