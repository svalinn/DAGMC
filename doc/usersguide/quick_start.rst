DAGMC Quick Start Guide
------

This quick start guide is intended to give a very brief concise flavor as to how to procede running DAGMC problems, it is assumed that you are familar with

* Cubit and the process of assigning volumes to groups, imprinting and merging
* That you have an installed MOAB,HDF5 and a DAGMC enabled MC code

If these things do not make sense, or you haven't read the remaining documentation, the `install guide <get_install.html>`_, `general workflow <workflow.html>`_ and
the `Unified Workflow <uw2.html>`_ are good places to start.

The general workflow is as follows;

1. Import your model into Cubit, assigning volumes to groups of the same material, i.e. adding all volumes that are to be steel to one group and so on
2. Simplify and defeature your CAD as your see fit
3. Imprint and merge the model, ```CRITICAL```
4. Export your model as an ACIS file, remember to check the export attributes box in the save window
5. Use dagmc_preproc, a tool found in the moab/bin directory, this is used to produce a faceted representation of the CAD model
6. Use make_watertight, a tool that can be downloaded from the Svalinn github site, to seal your geometry.
7. If you are using the UWUW workflow, use uwuw_preproc from the DAGMC/tools directory to mark up your model
8. Write the input deck for the code you want to run (you should've done this whilst you were faceting and sealing)
9. Run the problem as you see fit


Many models fail at the dagmc_preproc stage, due to problems with the CAD, overlapping volumes and surfaces, you must rectify these issues before continuing. Many downstream issues can be caused by lack of attention to detail at this stage. Lost particle problems can be fixed with make_watertight, we fully endorse the use of this tool.
