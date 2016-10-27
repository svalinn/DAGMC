The Cubit/Trelis plugin
=======================

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

Roadmap for the future
~~~~~~~~~~~~~~~~~~~~~~

TODO: this
