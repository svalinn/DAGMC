Installing the Cubit/Trelis Plugin
++++++++++++++++++++++++++++++++++

Linux Install Instructions
--------------------------
Obtain the `DAGMC Plug-in <>`_

   1.  unpack this tar file in the plugin directory of Trelis, e.g. /opt/Trelis-16.0/bin/plugins/
   2.  make a symbolic link in the plugin directory to the dagmc plugin itself. e.g. 
       ln -s dagmc/libdagmc_export_plugin.so .

OS/X Install Instructions
-------------------------
Obtain the `DAGMC Plug-in <>`_

   1.  unpack this tar file in the plugin directory of Trelis, e.g. /opt/Trelis-16.0/bin/plugins/
   2.  make a symbolic link in the plugin directory to the dagmc plugin itself. e.g. 
       ln -s dagmc/libdagmc_export_plugin.so .

Windows Install Instructions
--------------------------
Obtain the `DAGMC Plug-in <>`_

Using the Plug-in
--------------------
Load your geometry into Cubit/Trelis following the instructions in `geometry preparation<>`_, markup your 
geometry according to the workflow of your choice, for example `UWUW workflow <>`_, and Imprint and Merge. 
You can now export the geometry, ready for use in simulations using
::
   %> export dagmc <filename> [faceting_tolerance <faceting tolerance>] 
                              [length_tolerance <length tolerance>]
                              [normal_tolerance <normal tolerance>] 
			      [verbose] [fatal_on_curves]

Where faceting_tolerance, normal_tolerance, and length tolerance are optional arguments. For example, export
the currently loaded file to a file called "large_facets.h5m" with a faceting tolerance of 1.e-5 
::
   %> export dagmc "large_facets.h5m" faceting_tolerance 1.e-5 

