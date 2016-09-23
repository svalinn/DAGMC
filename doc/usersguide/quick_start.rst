Quick Start Guide
=================

This quick start guide is intended to give a very brief concise flavor as to how to proceed running
DAGMC problems.  It is assumed that you:

* are familiar with Trelis/Cubit and the process of assigning volumes to groups, imprinting and merging
  and how to proceed with one of the DAGMC workflows either UWUW for all the codes, or one of the specific
  MC code worklows
* have an installed MOAB, HDF5 and a DAGMC enabled MC code

If these things do not make sense, or you haven't read the remaining documentation, the `install guide <../install.html>`_
and `workflow <workflow/index.html>`_ are good places to start.

The general workflow is as follows;

  1. Import your model into Trelis/Cubit, assigning volumes to groups of the same material, e.g.
     adding all volumes that are to be steel to one group and so on
  2. Simplify and defeature your CAD as you see fit
  3. Imprint and merge the model
  4. Export your model to DAGMC
  5. Use make_watertight, a tool that is built alongside DAGMC, to seal your geometry
  6. If you are using the UWUW workflow, use uwuw_preproc from the DAGMC/bin directory to mark up your model
  7. Write the input deck for the code you want to run (Ideally you should've done this whilst you were faceting and sealing)
  8. Run the problem as you see fit

Many models fail at the "Export to DAGMC" stage due to problems with the CAD, such as overlapping volumes and surfaces.
You must rectify these issues before continuing, see `Using Trelis/Cubit for the DAGMC Workflow <workflow/cubit_trelis_workflow.html>`_ for guidance. Many downstream issues can be caused by lack of attention to detail
at this stage. Lost particle problems can sometimes be fixed with make_watertight; we fully endorse the use of this tool
, it is installed as part of the standard DAGMC tools.
