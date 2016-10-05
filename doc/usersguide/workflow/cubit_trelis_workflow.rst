Using Trelis/Cubit for the DAGMC Workflow
=========================================

The general workflow for the production of models for analysis using DAGMC
looks like that shown in the Figure below.

.. image:: general_workflow.png
   :height: 700
   :width:  600
   :alt: The general workflow for producing quality CAD for DAGMC models.

The general workflow for the production of DAGMC models is the following:

1. Ensure that units/sizes for DAGMC models are "cm"
2. Remove excessive detail (typically threads on bolts, combine washer stacks, etc.)
3. Inspect and resolve overlapping volumes
    - may need to scale model up to detect these
    - In Cubit:
        **validate vol all**
        **autoheal problem vols**
    - regularize may be needed (simplifies by removing splines, byproduct is reverses imprint)
    - Check and possibly remove small features (curves, surfaces, vols)
        - small areas (check "hydraulic" length and "regular" length)
            **group "smallsurfaces" add surface with area < 1.e-3**
        - small curves (check "hydraulic" length and "regular" length)
            **group "smallcurves" add curve with length <1.e-4**
        - small volumes
            volofvolcubit.py -- examine small volumes (vol<=0.0000 might need to autoheal)
4. create pre-imprint/merge table of volume of volumes (use volofvolcubit.py)
    - Imprint model, **imprint body all**
5. Merge model, in Cubit **merge all**
    - you may wish to add a merge tolerance, **merge tol 1e-6**
6. Validate model, in Cubit **validate vol all**
    - re-check for overlapping volumes and now also check for overlapping surfaces
    - create post imprint/merge table of volume of volumes (use volofvolcubit.py)
    - compare pre and post-imprint/merge model volume of volumes
    - inspect volumes in pre-imprint model with significant change in volume
    - re-check for small areas, curves, volumes
    - any problems in these steps repair volumes in pre-imprint model and go back to step 2
7. Facet model
    - **export dagmc "geom.h5m" faceting_tolerance 1.0e-4**
8. Seal model if possible (use make_watertight, does not always work)
9. Flood and/or transport particles in model
    - examine lost locations (use mklostvis.pl)
    - examine "leaks"/tunneling (can use a mesh tally to locate)
10. if lost particles or leaky repair the pre-imprint/merge model and go to step 2
11. if no lost or leaks, then transport is ok

Preparing Solid Models
~~~~~~~~~~~~~~~~~~~~~~

In theory, solid models can be prepared in any modeling software
system (e.g. SolidWorks, Pro/E, Catia, etc).  What is most important
about the choice of solid modeling system is the ability to export to
a format that can be imported by Trelis or Cubit, in particular:

    * ACIS (\*.sat)
    * STEP (\*.stp, \*.step, etc)

We have a strong preference towards ACIS due to its ability to retain 'imprint and
merge' information, and we have anecdotal evidence that ACIS files 
lead to a more successful model pipeline. Indeed we rely on Ansys Spaceclaim
to perform model cleaning and defeaturing, which can import many of the common
CAD formats and exports to ACIS.

There are a number of concepts to consider while preparing a solid
model; however, the most important considerations are small gaps and
overlaps that might exist in the model. These gaps and overlaps can
lead to rapid failure when running a DAGMC-based analysis. The
following steps are provided to help make a more robust model *before*
running your DAGMC-based analysis.

Be aware: obtaining a robust model may be an iterative and time
consuming process. In some cases, the validity of the model will
require running a DAGMC-based analysis and assessing whether or not
the model yielded expected results or a small enough number of lost
particles. If the results did not meet expectations, changes to the
model may be in order.

Knowing the model
-----------------

The first consideration to address is where the solid model originated
and for what purpose. In many instances, models constructed for
manufacturing purposes will have tolerances that are undesirable for
particle transport applications. For example, a gap might exist
between fuel pellets and the cladding wall for a PWR fuel rod. While
this is perfectly acceptable for an individual manufacturing the rod,
the gap could potentially present problems in a DAGMC-based
analysis, depending on how it is modeled.

Knowing who created the model and to what purpose provides a starting
point for preparing the model. If it was made with particle transport
in mind, then very little work may be needed; but as with the example
above, some models may require changes to accommodate the needs of a
DAGMC-based analysis.

Identifying weaknesses in the model
-----------------------------------

When assessing a model that is to be used for particle transport two
primary concerns must addressed. These concerns are:

    * Gaps
    * Overlaps

Gaps occur when the surfaces of two volumes/parts that should be in
contact are set apart from each instead of having coincident
surfaces. The size of the gap is generally unimportant, for most solid
modeling programs, a gap is a gap. The desired result is to have all
surfaces of volumes/parts to be coincident. If coincidence is not
achieved, particles may become lost when entering the region between
the surfaces.

Overlaps are found where two or more volumes/parts encroach upon the
same space. As with gaps, the magnitude of the overlapping volume is
usually unimportant.  When a particle enters a region of overlap, it
may not correctly determine which volume/part it resides in. If this
occurs, the particle may become lost.

Identifying gaps and overlaps may be difficult and time consuming;
however, some 3D modeling programs like SolidWorks have built in tools
to identify these occurrences. Rely on the modeling program to
identify these errors (the gaps and overlaps) and use the steps in the
next section to change, reduce and remove their effect on the model.

Modifying your model
--------------------

Once the gaps and overlaps in the model have been identified, the
three following methods may be used to change, reduce and remove their
effect on the model.

* Create "voided" geometries
* Modify volume/part dimensions
* Remove superfluous details

Each method is discussed in detail below:

As with the fuel rod example mentioned above, some geometries that are
'gaps' are also important. Instead of removing the gap entirely (by
changing the dimensions of the cladding or the fuel to force
coincidence), a new volume/part could be modeled that coincided with
the outer diameter of the fuel AND the inner diameter of the
cladding. Now a "voided" geometry occupies the previously unaccounted
for region. By specifying these "voided" geometries in a DAGMC-based
analysis, the physical importance of the region can be retained while
accomodating the requirement of having coincident surfaces.

Another method to resolve gaps and overlaps is to simply change the
dimensions of the volume/part (eg: making a dimension several cm
bigger or smaller to ensure coincidence surfaces). In many instances
this method could compromise the physics of the solution and is then
undesirable. However, in other instances, this solution is very
logical. One particularly significant example is if different volumes
were modeled with different unit systems. For example, one volume/part
might have been model in [in] while its neighbor was modeled in [cm];
while the surfaces may be nearly coincidence, rounding errors might
prevent coincidence from occurring. A simple change to one dimension
may hardly change the volume/part's characteristics yet result in
coincidence.

Finally, superfluous details may prevent a volume/part from coinciding
with its neighbors properly. A potential solution is to simply remove
the superfluous detail to simplfy the model and ensure the desired
surfaces are coincident. Some volumes/parts will inherently hurt the
model's effectiveness either due to its complex features or small
dimensions. A volume/part's effect on the model cannot truly be
assessed until a DAGMC-based analysis is run. This final method is
usually implemented in an attempt to reduce the number of lost particles
while maintaining the most important characteristics of the system.

*Note: Of all steps, the removal of superfluous details is the most
subjective and heavily dependent on the model's intended
application.*

Assessing your model
--------------------

Lost particles are undesirable; lost particles usually indicate
weaknesses and failures within the geometry. While the goal of the
DAGMC project is to guarantee that there will never be lost particles,
they can occur even on robust geometries.  It is up to the
user/analyst to determine what lost particle rate they consider
acceptable.  The UW-Madison group usually considers lost particle
rates that are less than 1/5,000,000 to be a threshold for most problems.
It is important to understand whether particles are being lost from an
important region of your phase space.

The implicit compliment is automatically generated by DAGMC upon loading a geometry;
it is composed of all the space that is not defined by the CAD geometry. It is often
convenient to not define all space in a given model, for example the space inside a
tokamak which is occupied by air or vacuum, or the water volume in a reactor. The
power of the implicit compliment lies in the fact that it is not a true CAD body
since it was never defined, but automatically defines all undefined space in the model.

Pre-processing Solid Models using Cubit/Trelis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Note: For large models, the steps described below can be very tedious
and time consuming.  To accelerate the process, an automated approach
is available for models that have been properly prepared in the native
solid modeling software.  This AutomatedCubitConversion process is
described elsewhere, but reading the information below will provided
the knowledge-base needed to understand the automation process.*

This section focuses on steps that are independent of the MC code used
for analysis. Additional steps for `DAG-MCNP5 <uw2.html>`_ and
`DAG-Tripoli4 <dag-tripoli4.html>`_ may be based on the instructions given here,
but are provided in the respective document links.

Importing the Solid Model
-------------------------

The first step in Cubit/Trelis is to import the generated solid
model. Depending on the complexity of the model, this step can take
several seconds up to a half an hour. As an initial user, it is
recommend to start with simple models and geometries to obtain a
better understanding of Cubit/Trelis.

Imprint and Merge
-----------------

For a DAGMC based analysis to work optimally, all of the surfaces must
be imprinted and merged. Imprinting creates a common surface
interface between touching volumes.  Merging then takes the touching
surfaces and makes them into one surface. The imprint operation is shown
below

.. image:: imprint_operation.png
   :height: 200
   :width:  600
   :alt: Imprint operations, results in the creation of additional surfaces.

To imprint, issue the following command:
::

    CUBIT> imprint body all

Should the imprint be successful, then the next step is to merge the
geometry. A schematic of what the merge operation achieves is shown
below,

.. image:: merge_operation.png
   :height: 250
   :width:  600
   :alt: Merge operations, results in the removal of identical surfaces.

Sometimes it may be important to specify a merge tolerance.
To set the tolerance and merge, issue the following commands:
::

    CUBIT> merge tol 5e-7
    CUBIT> merge all

This process can be very time consuming. For large models of several
thousand volumes, the imprint and merge steps can take several hours.
However, for small geometries (on the order of 100 volumes) the
process is rather quick.

.. _grouping-basics:

Grouping Volumes and Surfaces
-----------------------------

A DAGMC-based analysis allows a number of attributes of the geometry
to be defined within the geometry file. These characteristics
generally relate to the physical behavior of the volume, for example
its material definition or boundary conditions.

Before the discussion of specific attributes, the practice of
"grouping" needs to be explained. A group is essentially a collection
of volumes or surfaces that share a common attribute; the practical
usage of "grouping" will be explained in the next section.

The general format for creating/adding volumes to a group is:
::

    CUBIT> group "group.name" add vol/surf ...

For example, to create a group called "moderator" containing volumes
5, 6, 7, and 10, the following command would be used:
::

    CUBIT> group "moderator" add vol 5 to 8 10

Another example, shows that groups don't have to just contain
volumes, but can contain surfaces too. Below the group
"shield.boundary" is created with surfaces 16 and 37:
::

    CUBIT> group "shield.boundary" add surf 16 37

Due to the importance of using the ``group`` command reading the CUBIT
manual section on its full usage is highly recommended.

Production of the DAGMC Geometry
--------------------------------

Now that the geometry is ready for DAGMC we must export it. Using the
Cubit/Trelis plugin make this very straightforward, assuming that the
user has proceeded through the previous steps then all one must do is
use the export command, for example to produce a file called, geometry.h5m
with faceting tolerances and length tolerances of 1.0e-4 cm and 5.0 cm respectively
::

    CUBIT> export dagmc geometry.h5m faceting_tolerance 1.e-4 length_tolerance 5.0

The time taken to perform this step depends upon the complexity of the model, it could 
take seconds for very simple models to hours for very complex models. It is also possible
that faceting artifacts or failures could occur at this point, so monitor the output
of this command in the Cubit/Trelis command line. If issues due occurs, these should be addressed 
following the workflow listed above.

Finishing Up and Final Notes
----------------------------
Having prepared your model to completion with the appropriate groups created
, you can choose to save your model in various formats. Previously 
we recommended ACIS \*.sat files, but any format that reliably retains
imprortant metadata.  Recommended storage formats are ACIS, \*.Trelis or 
\*.cub files.

One should also use the `make_watertight <watertightness.html>`_ tool on the 
produced DAGMC \*.h5m file in order to completely seal your geometry, this 
should help prevent tolerance issues due to faceting.
