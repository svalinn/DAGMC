Installing the Cubit/Trelis Plugin
==================================

There is a common location that stores all the plugins that are available for use. Go
to the `DAGMC Plugins <http://go.wisc.edu/dagmc-trelis>`_ page and download the plugin appropriate
for your operating system.

Linux Install Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Obtain the Linux plugin from the link above. The following operations should be
performed as the root user (sudo).

1. Unpack tar ball in /opt/Trelis-16.0/bin/plugins
2. Change directory to /opt/Trelis-16.0/bin/plugins/dagmc
3. Run the install.sh script: ./install.sh

OS/X Install Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

Obtain the OS/X plugin from the link above. The following operations should be
performed as the root user (sudo).

1. If it does not already exist, create a folder called "plugins" in
   /Applications/Trelis-16.0/Contents/MacOS/
2. Unpack tar ball in /Applications/Trelis-16.0/Contents/MacOS/plugins
3. Change directory to /Applications/Trelis-16.0/Contents/MacOS/plugins/dagmc
4. Run the install.sh script: ./install.sh

Windows Install Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A plugin for Windows is currently being developed; it will be posted to the
above link when complete.

Using the Plugin
~~~~~~~~~~~~~~~~

Load your geometry into Cubit/Trelis following the instructions in the
`geometry preparation guide <workflow/cubit_trelis_workflow.html>`_, then markup your
geometry according to the workflow of your choice; for example,
`the UWUW workflow <workflow/uw2.html>`_, and then imprint and merge. You can then export the
geometry to a useable state for simulations with
::

   CUBIT> export dagmc <filename> [faceting_tolerance <faceting tolerance>]
                                  [length_tolerance <length tolerance>]
                                  [normal_tolerance <normal tolerance>]
                                  [verbose] [fatal_on_curves]

Where faceting_tolerance, normal_tolerance, and length tolerance are optional arguments. For example, to export
the currently loaded file to a file called "large_facets.h5m" with a faceting tolerance of 1.e-5, use
::

   CUBIT> export dagmc "large_facets.h5m" faceting_tolerance 1.e-5
