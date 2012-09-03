DAGMC User's Guide
========

This document is intended for users of the DAGMC toolkit and includes
instructions for installation and use with a variety of Monte Carlo
codes.

----

Getting & Installing the DAGMC Toolkit
----------------------------------------

DAGMC is a toolkit that provides direct geometry support to any Monte
Carlo radiation transport code.  The primary code used for development
and testing of DAGMC is `MCNP5 <http://mcnp-green.lanl.gov/>`_,
developed by `Los Alamos National Laboratory <http://www.lanl.gov>`_
and distributed by the `Radiation Safety Information Computing Center
<http://rsicc.ornl.gov>`_.  There has also been experience with MCNPX
(LANL) and Tripoli4 (CEA/Saclay).

These instructions describe the basic steps for downloading and
installing the software libraries that provide the DAGMC toolkit for
integration with Monte Carlo codes.  After this, code-specific
instructions will be give for each code. For access to the DAGMC
toolkit alone, the following combination of software
packages/libraries is necessary:

* `MOAB/DAGMC <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_
  (SVN: https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk)
    * `HDF5 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_
    * `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_ (SVN:
      https://svn.mcs.anl.gov/repos/ITAPS/cgm/trunk)
        * ACIS v16, or `CUBIT <http://cubit.sandia.gov>`_ v10.2

Installing the DAGMC Toolkit
--------

The following 4 steps are required to install the MOAB library,
including the DAGMC toolkit, for use in Monte Carlo radiation
transport tools.

1. Install `CUBIT <http://cubit.sandia.gov>`_ v10.2
2. Install `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_, using the --with-cubit options
3. Install `HDF5 <http://www.hdfgroup.org/HDF5/>`_
4. Install `MOAB <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_,
   using the --with-cgm and --with-hdf5 options (--with-netcdf may
   also be useful but not necessary)

Here are some assumptions/conventions that are used in these instructions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``
* all tarballs reside in user's home directory

If these do not apply to you, please modify your steps accordingly.

     *(For a shortcut to installing DAG-MCNP5.1.51 you may be able to use the DagmcBuildPackage .)*

Installing CGM
--------
::
    prompt%> mkdir -p $HOME/dagmc_bld/CGM/bld
    prompt%> cd $HOME/dagmc_bld/CGM

If installing from SVN repository:
::
    prompt%> svn co https://svn.mcs.anl.gov/repos/ITAPS/cgm/trunk
    prompt%> cd trunk
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s trunk src

If installing from a tarball, ``CGM-10.2.2.tar.gz``:
::
    prompt%> tar xzf ~/CGM-10.2.2.tar.gz
    prompt%> ln -s CGM-10.2.2 src

In all cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize --disable-debug \
              --with-cubit=/path/to/cubit  \
              --prefix=$HOME/dagmc_bld/CGM
    prompt%> make
    prompt%> make install

 1.A.ii Installing HDF5
-----------------------

::
    prompt%> mkdir $HOME/dagmc_bld/HDF5/bld
    prompt%> cd $HOME/dagmc_bld/HDF5
    prompt%> tar xzf ~/hdf5-1.8.3.tar.gz
    prompt%> ln -s hdf5-1.8.3 src
    prompt%> cd bld
    prompt%> ../src/configure --prefix=$HOME/dagmc_bld/HDF5
    prompt%> make
    prompt%> make install


.. _S1Aiii

 1.A.iii Installing MOAB
--------
::
prompt%> mkdir -p $HOME/dagmc_bld/MOAB/bld
prompt%> cd $HOME/dagmc_bld/MOAB


If installing from SVN repository:
::
prompt%> svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk
prompt%> cd trunk
prompt%> autoreconf -fi
prompt%> cd ..
prompt%> ln -s trunk src


If installing from a tarball, ``MOAB-3.99.tar.gz``:
::
prompt%> tar xzf ~/MOAB-3.99.tar.gz
prompt%> ln -s MOAB-3.99 src


In all cases:
::
prompt%> cd bld
prompt%> ../src/configure --enable-optimize --disable-debug \
          --with-cgm=$HOME/dagmc_bld/CGMA  \
          --with-hdf5=$HOME/dagmc_bld/HDF5 \
          --prefix=$HOME/dagmc_bld/MOAB
prompt%> make
prompt%> make install


----

.. _S1B

 1.B. Applying DAGMC to Specific Monte Carlo Codes
--------
.. _S1Bi

 1.B.i Installing DAG-MCNP5
--------
If you would like to use DAGMC with MCNP5, known as DAG-MCNP5, you will also need:
* MCNP5.1.51 source code from `RSICC <http://rsicc.ornl.gov>`_
* DAG-MCNP5.1.51 patch file from the UW-Madison

.. _S1Bia

 1.B.i.a Automatic Installation
--------
A package has been prepared that includes many of the requires software libraries and an automated build script.  Because the DAGMC team is not authorized to distribute `CUBIT <http://cubit.sandia.gov>`_ nor  `MCNP5.1.51 source code <http://mcnp.lanl.gov>`_, you must acquire those through the appropriate channels on your own.

Once you have both of those things, you should be able to use the DagmcBuildPackage to create a working install of DAG-MCNP5.1.51.

.. _S1Bib

 1.B.i.b Manual Installation
--------
The following steps are required to install DAG-MCNP5.  Most of these steps are described in more detail below.

1. Install the DAGMC Toolkit as described above
2. Download a copy of the patch file for your version of MCNP:
    * `MCNP5 v1.51 </software/dagmc.patch.5.1.51>`_
    * `MCNP4 v1.40 </software/dagmc.patch.5.1.40>`_
3. Apply the patch your copy of the MCNP5.1.51 source code
4. Build & install the patched version of MCNP5

Some assumptions/conventions:
* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``
* all tarballs reside in user's home directory
* MCNP5 source code is available in location ``$HOME/dagmc_bld/MCNP5``

**Apply DAGMC Patch to MCNP5 v1.51**

::
prompt%> cd $HOME/dagmc_bld/MCNP5/Source
prompt%> patch -p 0 < ~/dagmc_install.patch


**Build DAG-MCNP5 from modified code**

One of the easiest ways to build DAG-MCNP5 is directly using the ``makefile`` from the command-line.  To do this, you must know the ``makefile`` options to build a version of MCNP5 without DAGMC, usually in the form:

``prompt%> make build CONFIG``"seq plot gfortran" FC=gfortran MARCH=M64=

or similar.  Starting from these options, you can build DAG-MCNP5 from a patched source code with:

::
prompt%> make build CONFIG="seq plot gfortran dagmc" FC=gfortran MARCH=M64 \
             MOAB_DIR=$HOME/dagmc_bld/MOAB CUBIT_DIR=/path/to/cubit


If you are less familiar with building MCNP5 from the ``makefile`` you may want to use the interactive ``install`` script provided by LANL:

``prompt%> ./install``

Within the ``install`` program you will need to set the DAGMC build options:
* turn on DAGMC mode
* provide the path to MOAB: ``$HOME/dagmc_bld/MOAB``
* provide the path to CUBIT: ``/path/to/cubit``
From the main menu, choose ``C`` to compile.

Your executable should be available as ``$HOME/dagmc_bld/MCNP5/Source/src/mcnp5``.

.. _S1Bii

 1.B.ii Access to DAG-Tripoli4
--------
Tripoli4 is distributed by CEA/Saclay as a binary executable.  For access to DAG-Tripoli4, please contact `Jean-Christophe Trama <mailto:jean-christophe.trama@cea.fr>`_.

----
----

.. _S2

 2. Using DAGMC - The DAGMC Workflow
--------
The following steps are necessary to perform an analysis with most DAGMC versions of Monte Carlo radiation transport codes. (see below for more details on each step)
1. create and/or prepare a solid model using the CAD/solid-modeling tool of your choice
2. pre-processing of the solid model geometry in CUBIT
    * assign materials and densities
    * define boundary conditions
    * imprint & merge geometry
    * export that model in the ACIS format
3.  prepare native MC code input file
    * understanding additional ``dagmc`` parameters
4. Run MC code, possibly with additional command line parameters

----

.. _S2A

 2.A. Preparing Solid Models
--------
In theory, solid models can be prepared in any modeling software system (e.g. SolidWorks, Pro/E, Catia, etc).  What is most important about the choice of solid modeling system is the ability to export to a format that can be imported by CUBIT, in particular:
* ACIS (*.sat)
* STEP (*.stp, *.STEP, etc)

There are a number of concepts to consider while preparing a solid model; however, the most important considerations are small gaps and overlaps that might exist in the model. These gaps and overlaps can lead to rapid failure when running a DAGMC-based analysis. The following steps are provided to help make a more robust model *before* running your DAGMC-based analysis.

Be aware: obtaining a robust model may be an iterative and time consuming process. In some cases, the validity of the model will require running a DAGMC-based analysis and assessing whether or not the model yielded expected results or a small enough number of lost particles. If the results did not meet expectations, changes to the model may be in order.

.. _S2Ai

 2.A.i: Knowing the model
--------
The first consideration to address is where the solid model originated and for what purpose. In many instances, models constructed for manufacturing purposes will have tolerances that are undesirable for particle transport applications. For example, a gap might exist between fuel pellets and the cladding wall for a PWR fuel rod. While this is perfectly acceptable for an individual manufacturing the rod, the gap could potentially cause present problems in a DAGMC-based analysis, depending on how it is modeled.

Knowing who created the model and to what purpose provides a starting point for preparing the model. If it was made with particle transport in mind, then very little work may be needed; but as with the example above, some models may require changes to accommodate the needs of a DAGMC-based analysis.

.. _S2Aii

 2.A.ii: Identifying weaknesses in the model
--------
When assessing a model that is to be used for particle transport two primary concerns must addressed. These concerns are:

    * Gaps
    * Overlaps

Gaps occur when the surfaces of two volumes/parts that should be in contact are set apart from each instead of having coincident surfaces. The size of the gap is generally unimportant, for most solid modeling programs, a gap is a gap. The desired result is to have all surfaces of volumes/parts to be coincident. If coincidence is not achieved, particles may become lost when entering the region between the surfaces.

Overlaps are found where two or more volumes/parts encroach upon the same space. As with gaps, the magnitude of the overlapping volume is usually unimportant.   When a particle enters a region of overlap, it may not correctly determine which volume/part it resides in. If this occurs, the particle may become lost.

Identifying gaps and overlaps may be difficult and time consuming; however, some 3D modeling programs like SolidWorks have built in tools to identify these occurrences. Rely on the modeling program to identify these errors (the gaps and overlaps) and use the steps in the next section to change, reduce and remove their effect on the model.

.. _S2Aiii

 2.A.iii: Modifying your model
--------
Once the gaps and overlaps in the model have been identified, the three following methods may be used to change, reduce and remove their effect on the model.

* Create "voided" geometries
* Modify volume/part dimensions
* Remove superfluous details

Each method is discussed in detail below:

As with the fuel rod example mentioned above, some geometries that are 'gaps' are also important. Instead of removing the gap entirely (by changing the dimensions of the cladding or the fuel to force coincidence), a new volume/part could be modeled that coincided with the outer diameter of the fuel AND the inner diameter of the cladding. Now a "voided" geometry occupies the previously unaccounted for region. By specifying these "voided" geometries in a DAGMC-based analysis, the physical importance of the region can be retained while accomodating the requirement of having coincident surfaces.

Another method to resolve gaps and overlaps is to simply change the dimensions of the volume/part (eg: making a dimension several cm bigger or smaller to ensure coincidence surfaces). In many instances this method could compromise the physics of the solution and is then undesirable. However, in other instances, this solution is very logical. One particularly significant example is if different volumes were modeled with different unit system. For example, one volume/part might have been model in [in] while its neighbor was modeled in [cm]; while the surfaces may be nearly coincidence, rounding errors might prevent coincidence from occurring. A simple change to one dimension may hardly change the volume/part's characteristics yet result in coincidence.

Finally, superfluous details may prevent a volume/part from coinciding with its neighbors properly. A potential solution is to simply remove the superfluous detail to simplfy the model and ensure the desired surfaces are coincident. Some volumes/parts will inherently hurt the model's effectiveness either due to its complex features or small dimensions. A volume/part's effect on the model cannot truly be assessed until a DAGMC-based analysis is run. This final method is usually implemented in attempt to reduce the number of lost particles while maintaining the most important characteristics of the system.

*Note: Of all steps, the removal of superfluous details is the most subjective and heavily dependent on the model's intended application._

.. _S2Aiv

 2.A.iv: Assessing your model
--------
Lost particles are undesirable; lost particles usually indicate weaknesses and failures within the geometry. While the goal of the DAGMC project is to guarantee that there will never be lost particles, they can occur even on robust geometries.  It is up to the user/analyst to determine what lost particle rate they consider acceptable.  The UW-Madison group usually considers lost particle rates that are less than 1/50,000 to be a threshold for most problems.  It is important to understand whether particles are being lost from an important region of your phase space.

[Insert note on the implicit complement here]

----

.. _S2B

 2.B. Pre-processing Solid Models using CUBIT
--------
_Note: For large models, the steps described below can be very tedious and time consuming.  To accelerate the process, an automated approach is available for models that have been properly prepared in the native solid modeling software.  This AutomatedCubitConversion process is described elsewhere, but reading the information below will provided the knowledge-base needed to understand the automation process.*

This section focuses on steps that are independent of the MC code used for analysis.  Additional steps for `DAG-MCNP5 <#S2Di>`_ and `DAG-Tripoli4 <#S2Dii>`_ may be based on the instructions given here, but are provided in separate parts of section 2.D below.

.. _S2Bi

 2.B.i Importing the Solid Model
--------
The first step in CUBIT is to import the generated solid model. Depending on the complexity of the model, this step can take several seconds up to a half an hour. As an initial user, it is recommend to start with simple models and geometries to obtain a better understanding of CUBIT.

.. _S2Bii

 2.B.ii Imprint and Merge
--------
For a DAGMC-based analysis to work properly, all of the surfaces must be imprinted and merged.  Imprinting creates a common surface interface between touching volumes.  Merging then takes the touching surfaces and makes them into one surface.

To imprint, issue the following command

  imprint body all

Should the imprint be successful, then the next step is to merge the geometry. Sometimes it may be important to specify a merge tolerance.  To set the tolerance and merge, issue the following commands:

  merge tol 5e-7
  merge all

This process can be very time consuming. For large models of several thousand volumes, the imprint and merge steps can take up to three hours. However, small geometries (on the order of 100 volumes) the process is rather quick.


.. _S2Biii

 2.B.iii Grouping Volumes and Surfaces
--------
A DAGMC-based analysis allows a number of attributes of the geometry to be defined within the geometry file. These characteristics generally relate to the physical behavior of the volume, for example its material definition or boundary conditions.

Before the discussion of specific attributes, the practice of "grouping" needs to be explained. A group is essentially a collection of volumes or surfaces that share a common attribute; the practical usage of "grouping" will be explained in the next section.

The general format for creating/adding volumes to a group is:

  group "group.name" add vol/surf ...

For example, to create a group called "moderator" containing volumes 5, 6, 7, and 10, the following command would be used:

  group "moderator" add vol 5 to 8 10

Another example, shows that groups don't have to be just contain volumes, but can contain surfaces too. Below the group "shield.boundary" is created with surfaces 16 and 37:

  group "shield.boundary" add surf 16 37

Due to the importance of using the ``group`` command reading the CUBIT manual section on its full usage is highly recommended.


.. _S2Biv

 2.B.iv Finishing Up and Final Notes
--------
Before exporting, it is vital to set attributes on.  This saves the absolute volume and surface IDs as well as any group specifications.  Failing to do this will result in fatal errors.  Make sure to type the following:

  set attribute on

Finally, export the file as an ACIS file with a .sat extension.  If you are using a version of CUBIT newer that v10.x, be sure to set the ACIS geometry level back to version 16:

  set geom version 1600

For the remainder of this documentation, the geometry file will be referred to as "geom.sat". Also, as noted before, the CUBIT conversion process can be automated as described on the follow webpage: AutomatedCubitConversion.

----

.. _S2C

 2.C. Additional Input Parameters for DAGMC
--------
DAGMC introduces a number of new input variables that control the behavior of DAGMC when it interacts with your geometry.  This section describes the conceptual behavior of those parameters and later sections will describe their specific implementation and syntax within each code.

.. _S2Ci

 2.C.i Geometry File (basic)
--------
**required**

**Default: none**

This file contains the geometry that has been created and pre-processed using the workflow described above.  This file can be either an ACIS geometry (usually with a ``.sat`` file extension) or a MOAB facet file (usually with a ``.h5m`` file extension).

.. _S2Cii

 2.C.ii Faceting Tolerance (basic)
--------
**optional**

**Default: 0.001**

One of the first steps performed by DAGMC is to generate a faceted representation of your solid model.  To ensure a faithful representation of the model, the facets are constrained such that all points on each facet are within a tolerance of the nearest points on the exact surface representation.  A smaller tolerance results in a more faithful representation of the surface at the penalty of producing more facets.  The user can control the faceting tolerance using when they invoke their simulation, either on the command line or in the input file, depending on the MC code being used for the analysis.  This option only has an effect with the geometry file is a solid model and not when it is a facet file.

.. _S2Ciii

 2.C.iii Facet File (basic)
--------
**optional**

**Default: none**

For some models, the initial processing can be time consuming.  When reading a solid model geometry, this option causes a processed file to be written that can be used on subsequent analyses.  This file will be processed with the facet tolerance as defined above.  This facet tolerance cannot be changed when the file is reused.

.. _S2Civ

 2.C.iv Overlap Thickness (advanced)
--------
**optional**

**Default: 0.0**

Often CAD geometry has small overlaps between adjacent volumes due to file translation or imprecise modeling. The particle tracking algorithm can accommodate small overlaps if the maximum thickness of overlaps is approximately known.

.. _S2Cv

 2.C.v Source Cell Treatment (intermediate)
--------
**optional**

**Default: on (same behavior as native code)**

The implementation of this option is specific to the Monte Carlo code being used.  Please refer to the documentation for your underlying Monte Carlo code.

.. _S2Cvi

 2.C.vi Use CAD geometry (advanced)
--------
**optional**

**Default: off**

When this option is turned on, the ray-firing process finds the intersection with the CAD-based solid model itself, and not just with the faceted representation of that model.  The facet-based ray-firing solution is used as an initial guess to solve for the point on the actual CAD surface where the ray-surface intersection takes place.  This option is only available when the DAGMC toolkit has been linked to the ACIS geometry libraries directly and not when it has been linked via CUBIT.

.. _S2Cvii

 2.C.vii Use Distance Limit (experimental)
--------
**optional**

**Default: off**

This option allows a previously determined distance to the next collision to be used to accelerate the search for ray-surface intersections.  Any candidate ray-surface intersection that is at a distance beyond the known distance to collision will be rejected, including bounding box tests in the OBB tree search.

----

.. _S2D

 2.D. Monte Carlo Code-Specific Steps
--------
There are three varieties of code-specific steps:
1. defining attributes of the geometry using groups in CUBIT
2. defining DAGMC runtime parameters using input file syntax
3. changes to the command-line

In this section, these steps are described for each of the supported Monte Carlo codes.

.. _S2Di

 2.D.i. DAG-MCNP5
--------
.. _S2Dia

 2.D.i.a. Geometry Metadata
--------
In DAG-MCNP5, the geometry file can be used to define material and density assignments, boundary conditions, and tally locations.

    Assigning Materials & Densities
--------
The generic workflow description outlines the `principles of grouping <#S2Biii>`_, but a specific naming convention is required for DAG-MCNP5. To define materials, both the MCNP material ID and density must be provided in the group name. The format for the group name is as follows:

  mat_[matid]_rho_[density]

``[matid]`` is replaced by the material ID that will be specified in the MCNP input file.  ``[density]`` is replaced by either the atomic density or the mass density.  Like the MCNP cell cards, positive values are atomic densities in [atoms/barn-cm] and negative values are mass densities in [g/cc].

For example, suppose UO2 is material 7 in the problem, has an atomic density of 0.0223 and volumes 4 through 18 consist of this material.  To assign materials to these volumes, the following command would be used:

  group "mat_7_rho_0.0223" add vol 4 to 18

*Note: If a volume is not assigned to a specific group, when run in DAGMC it will be treated as a void; the material for that cell will be zero. This can actually become a fairly useful debugging tool to identify volumes that were not assigned to their appropriate group.*

If you would like to assign a material to the explicit complement, you can use the same mechanism, but add the string ``_comp`` to the end of the group name.  Since DAGMC only recognizes those groups that contain an entity, it is necessary to add a volume to this group, recognizing that this material will NOT be assigned to that volume.  (_It is often convenient to assign the graveyard volume (see below) to the implicit complement material group to minimize confusion._)  For example, if you would like the explicit complement to be modeled as material 9 with density 1 g/cc:

   create group "mat_9_rho_-1_comp"

    Defining Boundary Conditions
--------
There are two general classes of boundary condition supported by DAG-MCNP5. a vacuum boundary and reflecting surfaces, and they are implemented in different ways.

* **Defining the "graveyard": vacuum boundaries**

A typical usage of MCNP5 includes a volume that extends to infinity with an importance of 0 that bounds the active volumes of interest.  Since solid models cannot include any infinite volumes, it is necessary to place a finite volume of importance 0 to define the problem boundary. You will need to surround the entire geometry with a shell of finite thickness, known as the "graveyard".  Any geometric shape can be used for this; however, a cubic shell is preferred.  This shell volume will represent the outside world and will be the cell where all of the particles are terminated (thus it will have an importance of zero).

To create this volume create two volumes in CUBIT with the same shape, same center, and one slightly larger than the other.  Subtract the smaller from the larger.  The remaining volume is the graveyard.

Like the material definitions and boundary conditions discussed in the previous section. The graveyard is defined by assigning it a specific group name, one of the following keywords:

  graveyard
  outside.world
  rest.of.world

Consider a geometry with 99 volumes that all fit within a cube centered at the origin with side-length 99 cm.  To create a graveyard for this problem in CUBIT, you could issue the following commands:
::
cubit_prompt> create brick x 100
cubit_prompt> create brick x 105
cubit_prompt> subtract vol 100 from vol 101
cubit_prompt> group "graveyard" add vol 102


When DAG-MCNP5 is run, the importance of volume 102 (or any other volumes included in the group) will be set to zero. (_Note: this assumes that the two ``create brick`` commands generate volumes numbered 100 and 101, respectively, and that the Boolean subtraction results in a new volume number 102.

If you have boundary conditions (reflecting, white, or periodic) it is not required that you surround them with the bounding volume, but is not incorrect to do so.  Only areas where particles should escape need to be enclosed.  However, it is often easiest to simply create a single graveyard that covers all directions and volumes of the system.

* **Surface boundary conditions: reflection**

Surface boundary conditions are similarly enforced by specifying a group name. This type of attribute (surface boundary condition) is only required if reflective or white boundary conditions are used in the problem.  If not, this section may be skipped.  *Note that periodic boundary conditions are not yet supported.*

Specifying reflecting and white boundary conditions are fairly straightforward.  The group names for reflecting and white are respectively:

  spec.reflect
  white.reflect

Suppose surfaces 10 and 11 are reflecting boundary conditions.  To specify these as reflecting surfaces, the following group would be created:

  group "spec.reflect" add surf 10 11

*:Tally Assignments
--------
It is also possible, although not required, to specify tallies in the geometry.  The general form for adding this meta-data is to create a group of volumes or surfaces and encode the meta-data in the names of those groups.

The user has the option of specifying tallies in the geometry directly.  It is still possible to specify tallies in the MCNP input file, however, the user has to make sure that the tally indices are not duplicated lest a fatal error will occur.  Tallies are specified as group names in the following format:

  tally_[CUBIT tally ID]_[tally type keyword]_[particles]

The ``[CUBIT tally ID]`` field is an integer from 0 to 99.  Different tally types may have the same CUBIT ID and are still consistent.  The tally number in MCNP is 10 times the CUBIT ID plus the tally type index (e.g. 4 for cell flux tallies).

The ``[tally type keyword]`` is one of the following for each type of tally:

|Tally Type|tally type keyword|
|f1	|surf.current|
|f2	|surf.flux|
|f4	|cell.flux|
|f6	|cell.heating|
|f7	|cell.fission|
|f8	|pulse.height|

Also *tallies (the tally result times the incident particle energy) are possible by placing an "e" before the tally type.  So to make a *f2 tally, the keyword would be esurf_flux.  Pulse height (f8) tallies have the option to include charge as well.  This is done by placing a "q" before the keyword as in qpulse_height.

The ``[particles]`` tag is a string stating which particles will be tallied.  To tally both photons and neutrons, set the tag to "np".  The default is neutrons only.  Should this be tag be omitted, only neutrons will be tallied.

Some CUBIT commands to do tallies:

  group "tally_0_surf.current" add surf 1 to 4
  group "tally_0_cell.flux_p" add vol 7
  group "tally_1_ecell.heating_np" add vol 2 6
  group "tally_6_cell.heating_n" add vol 2 6
  group "tally_7_cell.flux_p" add vol 1 to 3
  group "tally_12_pulse.height_p" add vol 10 to 14
  group "tally_14_qpulse.height_p" add vol 10 to 14

The above are equivalent to following MCNP definitions:

  f1:n 1 2 3 4 T
  f4:p 7 T
  *f16:n,p 2 6 T
  f66:n 2 6 T
  f74:p 1 2 3 T
  f128:p 10 11 12 13 14 T
  +f148:p 10 11 12 13 14 T

*(Note: the current convention is to always add a tally bin for the total across all cells/volumes.)*

.. _S2Dib

 2.D.i.b. Preparing the DAG-MCNP5 Input File
--------
The DAG-MCNP5 input file contains only the data cards section of a standard MCNP5 input file.  There are no cell or surface cards included in the input file.

In addition to many other MCNP5 data cards, it is important to define the materials that have been assigned in step 2.D.i.a above and any tally modifiers, as desired, for the tallies defined in step 2.D.i.a above.

A new data card has been added to DAG-MCNP5 to define parameters for the DAGMC geometry capability.  These parameters are `described conceptually above <#S2C>`_.

::
Form: dagmc  keyword1=value   keyword2=value

      Keywords: check_src_cell usecad distlimit overlap_thickness

      overlap_thickness `` allows particle tracking through small overlaps
                          {real} [default``0.0]
      check_src_cell    `` behavior of CEL variable in SDEF card
                          on  [default] standard interpretation for CEL variable: source rejection
                          off           no cell rejection - assume that sampled position is in cell CEL
      usecad            `` toggle usage of solid model geometry
                          off [default] ray-tracing limited to facets
                          on            ray-tracing performed on solid model geometry surfaces
      distlimit         `` toggle usage of flight distance sampled from physics
                          to accelerate ray-tracing search
                          off [default] do not use physics flight distance
                          on            do use physics flight distance


.. _S2Dic

 2.D.i.c. Running DAG-MCNP5
--------
Running DAG-MCNP5 is identical to running the standard MCNP5, but a few new keywords have been added to the command-line to specify the necessary files.

``gcad=<geom_file>=
     (required) The ``geom_file`` is the geometry file that contains your geometric model, either in the ACIS (*.sat) format or the MOAB (*.h5m) format.  If this entry is not present, DAG-MCNP5 will assume that it is running in standard MCNP5 mode.  This runtime parameter is described in more detail above.

``ftol``<faceting_tolerance>=
     (optional) [default: 1e-3] This is a real number that provides guidance to the faceting engine regarding the maximum distance between a facet and the surface it is representing.  It is only used when reading an ACIS (*.sat) ``geom_file``.  When reading a MOAB (*.h5m) file, the facets have already been generated and this setting is ignored.  This runtime parameter is described in more detail above.

``fcad``<facet_file>
     (optional) The ``facet_file`` is written by DAG-MCNP5 in the MOAB (*.h5m) format.  When an ACIS file is read by DAG-MCNP5, a number of pre-processing and initialization steps are necessary.  Since these can be time consuming, the user has the option to create a ``facet_file`` the first time that they use a geometry and then use that ``facet_file`` with the ``gcad`` keyword in subsequent uses.  This runtime parameter is described in more detail above.


``lcad``<log_file>=
     (optional) The ``log_file`` is a skeleton of an MCNP file for the cells and surfaces in your geometry.  This file is created by DAG-MCNP5 to communicate the material assignments, boundary conditions, and tallies that you defined in your geometry.  If you give a name other than the default (=lcad=) for this file on the command-line, that file will be used instead of the one generated automatically by DAG-MCNP5.  This is useful to make small changes to your material assignments and/or importances, but **can not** be used to change the geometry.  It is up to the user to ensure that the ``log_file`` being used corresponds to the geometry file in question.  This runtime parameter is unique to the DAG-MCNP5 implementation of DAGMC.

.. _S2Dii

 2.D.ii. DAG-Tripoli4
--------
.. _S2Diia

 2.D.ii.a. Geometry Metadata
--------
The current version of DAG-Tripoli4 allows the definition of material compositions and boundary conditions in the geometry.

    Assigning Materials & Densities
--------
The generic workflow description outlines the `principles of grouping <#S2Biii>`_, but a specific naming convention is required for DAG-Tripoli4. To define materials, the Tripoli composition name must be provided in the group name. The format for the group name is as follows:

  comp_[compname]

``[compname]`` is replaced by the composition name that will be specified in the Tripoli4 input file.  ``[density]`` is replaced by either the atomic density or the mass density.  Like the MCNP cell cards, positive values are atomic densities in [atoms/barn-cm] and negative values are mass densities in [g/cc].

For example, suppose UO2 is composition ``Uoxide`` in the problem and volumes 4 through 18 consist of this material.  To assign materials to these volumes, the following command would be used:

  group "comp_Uoxide" add vol 4 to 18

*Note: If a volume is not assigned to a specific group, when run in DAGMC it will be treated as a void; the material for that cell will be zero. This can actually become a fairly useful debugging tool to identify volumes that were not assigned to their appropriate group.*

If you would like to assign a material to the explicit complement, you can use the same mechanism, but add the string ``_comp`` to the end of the group name.  Since DAGMC only recognizes those groups that contain an entity, it is necessary to add a volume to this group, recognizing that this material will NOT be assigned to that volume.  (_It is often convenient to assign the graveyard volume (see below) to the implicit complement material group to minimize confusion._)  For example, if you would like the explicit complement to be modeled with the composition named "air":

   create group "comp_air_comp"

    Defining Boundary Conditions
--------
There are two general classes of boundary condition supported by DAG-Tripoli4. a vacuum boundary and reflecting surfaces, and they are implemented in different ways.

* **Defining the "graveyard": vacuum boundaries**

A vacuum boundary condition is typically defined in Tripoli4 by simply having a surface with not defined volume on its other side.  Since DAGMC's implicit complement is also defined this way, it is not possible to use this convention in DAG-Tripoli4.  Instead, volumes are created on the vacuum side of those same surfaces, but those volumes are placed in a special group, known as the "graveyard" to change their behavior.  Any geometric shape can be used for this; however, a cubic shell is often preferred.  This shell volume will represent the outside world and will be the cell where all of the particles are terminated (thus it will have an importance of zero).

[need more on implicit complement in context of Tripoli]

To create this volume create two volumes in CUBIT with the same shape, same center, and one slightly larger than the other.  Subtract the smaller from the larger.  The remaining volume is the graveyard.

Like the material definitions and boundary conditions discussed in the previous section. The graveyard is defined by assigning it a specific group name, one of the following keywords:

  graveyard
  outside.world
  rest.of.world

Consider a geometry with 99 volumes that all fit within a cube centered at the origin with side-length 99 cm.  To create a graveyard for this problem in CUBIT, you could issue the following commands:
::
cubit_prompt> create brick x 100
cubit_prompt> create brick x 105
cubit_prompt> subtract vol 100 from vol 101
cubit_prompt> group "graveyard" add vol 102


When DAG-Tripoli4 is run, the surfaces of volume 102 (or any other volumes included in the group) will be defined as having only one volume - the one on the OTHER SIDE relative to voluem 102. (_Note: this assumes that the two ``create brick`` commands generate volumes numbered 100 and 101, respectively, and that the Boolean subtraction results in a new volume number 102.

If you have boundary conditions (reflecting, white, or periodic) it is not required that you surround them with the bounding volume, but is not incorrect to do so.  Only areas where particles should escape need to be enclosed.  However, it is often easiest to simply create a single graveyard that covers all directions and volumes of the system.

* **Surface boundary conditions: reflection**

Surface boundary conditions are similarly enforced by specifying a group name. This type of attribute (surface boundary condition) is only required if reflective or white boundary conditions are used in the problem.  If not, this section may be skipped.  *Note that periodic boundary conditions are not yet supported.*

Specifying reflecting and white boundary conditions are fairly straightforward.  The group names for reflecting and white are respectively:

  spec.reflect
  white.reflect

Suppose surfaces 10 and 11 are reflecting boundary conditions.  To specify these as reflecting surfaces, the following group would be created:

  group "spec.reflect" add surf 10 11

.. _S2Diib

 2.D.ii.b. DAGMC Runtime Parameters
--------
The DAGMC-Tripoli input file is formatted just like any other Tripoli input file but using the ``DAGMC_GEOMETRY`` block to indicate the geometry.  This block has the following parameters:

::
<geometry_filename>  This must be the first parameter
facet_tol <double faceting tolerance> (optional: default=0.001)
facet_file <string faceting filename> (optional)
check_src_cell <"off"|"false"|"no"> (optional: default=on)
usecad <"on"|"true"|"yes"> (optional: default=off)
distlimit <"on"|"true"|"yes"> (optional: default=off)
tolerance <double ray firing tolerance> (optional: default=1e-8)


These parameters are `defined in more detail above <#S2C>`_.  In addition to many other Tripoli input blocks, it is important to define the material compositions that have been assigned in step 2.D.ii.a.

.. _S2Diic

 2.D.ii.c. Running DAGMC-Tripoli
--------
Running DAGMC-Tripoli is identical to running the standard Tripoli.

----
----

.. _S3

 3. Frequently Asked Questions
--------
`Q1. How does the implicit complement work? <#Q1>`_

`Q2. Can I automate the CUBIT conversion process? <#Q2>`_

----

.. _Q1

 Q1:
--------*How does the implicit complement work?*

**A1** Since your geometry has been imprinted and merged, all surfaces fall into one of two categories:
1. surfaces that are shared by two volumes
2. surfaces that are used in only one volume

All surfaces in this second group are, by definition, on the boundary of the complement region since they have an explicitly defined region on one side and nothing on the other.  The implicit complement is formed automatically from this list of surfaces.

From the point of view of MCNP, and extra cell is created to represent the implicit complement.  Your output file will include one entry for each volume in your geometry, including the graveyard, plus one extra entry representing the complement.

.. _Q2

 Q2:
--------*The CUBIT conversion process is tedious and frustrating. Is there a way to avoid or automate most of this work?*

**A2** Yes, a script has been written that will automate all the necessary steps to perform the CUBIT conversion process with limited user input. Information about this script can be found on the AutomatedCubitConversion website.

----
----

.. _S4

 4. Related Tools
--------
* Automatic generation of graveyard
* Command-line conversion from ACIS (*.sat) to MOAB (*.h5m)
* Visualization of geometry in VisIt
* Visualization of mesh tallies in VisIt
* Converting MCNP geometry to CAD
* Automated CUBIT conversion script (AutomatedCubitConversion)

----
----
_[ HomePage - CnergSoftware - DirectAcceleratedGeometryMonteCarlo ]_
