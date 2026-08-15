Code-specific steps for FluDAG
==============================

There are several varieties of code-specific steps:

1.  Defining attributes of the geometry using groups in Cubit
2.  Producing material assignments in FLUKA input format from the h5m file, with
    the help of FluDAG
3.  Preparing the FLUKA input file for running with DAGMC
4.  Inserting the material assignments into the FLUKA input deck

Geometry metadata
~~~~~~~~~~~~~~~~~

In FluDAG, the geometry file can be used to define material assignments.
Eventually we would like to add the capability to define boundary conditions and
tally locations.

Assigning materials and densities
---------------------------------

The generic workflow description includes details on :ref:`grouping-basics`,
but a specific naming convention is required for FluDAG. To define materials,
the FLUKA material name must be provided in the group name. The group name must
be that of a valid FLUKA material; for example, any of the default materials or
a valid COMPOUND name. Note that you are limited to 8 characters, as in FLUKA.
The format for the group name is as follows:
::

    mat:[material_name]

For example, suppose we wish to add volumes 1 through 5 to a group that defines
the material to be iron. The following command would be used.
::

    CUBIT> group "mat:IRON" add volume 1 to 5

This will produce in the input file:
::

    ASSIGNMA        IRON         1
    ASSIGNMA        IRON         2
    ASSIGNMA        IRON         3
    ASSIGNMA        IRON         4
    ASSIGNMA        IRON         5

Compounds are also supported by FluDAG; for example, if we wish to have volume 6
belong to a group whose material name is STAINLESS then we can can use
::

    CUBIT> group "mat:STAINLESS" add volume 6

This will produce in the input file:
::

    MATERIAL                                        26                    STAINLES
    *...+....1....+....2....+....3....+....4....+....5....+....6....+....7...
    ASSIGNMA    STAINLES         6

Some notes:
    * Material names longer than 8 characters are truncated to the first 8
      characters.
    * There are several predefined material names in FLUKA, and they are
      appropriately treated by FluDAG.
    * All volumes must belong to a group. If any are not, FluDAG will not assign
      material information for them.

The implicit complement is automatically assigned the value 1 + the id of the
highest-numbered volume. You can easily modify what material property that you
would like the implicit compliment (or any other volume) to have by changing the
material on the ASSIGNMAT card.

Defining the graveyard
----------------------

Typical FLUKA models contain a zero-importance volume that bounds the volumes of
interest and extends to infinity. Since solid models cannot include infinite
volumes, it is necessary to define a finite zero-importance volume around your
model to define the problem boundary. This is done by surrounding the entire
geometry with a shell of finite thickness; this is known as the "graveyard." Any
geometric shape can be used for this, but a cubic shell is preferred. The
graveyard represents the outside world, and any particle that enters it will be
terminated.

To create a graveyard volume, create two volumes in Cubit with the same shape
and same center with one slightly larger than the other, making sure that both
bound the entire problem geometry. Then, subtract the smaller one from the
larger one. The remaining volume is the graveyard.

To indicate to FLUKA that a given volume is the graveyard volume, you must
assign it to a group with:
::

    CUBIT> group "mat:BLCKHOLE" add volume X

For example, consider a geometry containing 99 volumes, all of which fit inside
a cube of side length 99 cm centered at the origin. The following commands would
create a valid graveyard for this problem:
::

    CUBIT> create brick x 100             # This will be volume 100
    CUBIT> create brick x 105             # This will be volume 101
    CUBIT> subtract vol 100 from vol 101  # This will produce volume 102
    CUBIT> group "mat:BLCKHOLE" add vol 102

When FluDAG is run, the importance of any graveyard volumes will be set to zero.

Scoring assignments
-------------------

DAGMC does not currently support scoring assignments through group names. Users
must add these to the FLUKA input deck manually.

Preparing the FluDAG input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FluDAG input file will look almost identical to the originating FLUKA input
file. The exception will be the removal of all data between the cards
``GEOBEGIN`` and ``GEOEND``, i.e. all native FLUKA input data. The last entry
on the line of ``GEOBEGIN`` should be ``FLUGG``.

For example, the most simple valid FLUKA geometry is as follows:
::

    GEOBEGIN                                                              COMBNAME
        0    0
    SPH S1         0.0 0.0 0.0 50.0
    CELL1        5 +S1
    CELL2        5 -S1
    GEOEND

To run this geometry with FluDAG, remove all data between ``GEOBEGIN`` and
``GEOEND``, and switch the last entry to ``FLUGG``:
::

    GEOBEGIN                                                              FLUGG
    GEOEND

Running FluDAG
~~~~~~~~~~~~~~

Running FluDAG bears some similarity to running ``FLUGG``. The first step is to
create the CAD geometry of the problem you wish to run. In order to produce the
material assignment data from the CAD geometry we must first
:ref:`facet the file <geom_production>` using the Cubit plugin. Using the
subsequently-defined geometry file, the user must produce the ``mat.inp`` file.
::

    $ /path/to/fludag/executable/mainfludag --dagmc geom.h5m

This will load the named .h5m file and produce the material assignment
information. This information should then be pasted into the FLUKA input file,
and any adjustments that need to be made should be made; for example, adding the
density of non standard materials, or adding your scoring information. **Please
note that users must always include the additional material and compound
information themselves and take responsibility to ensure that the FLUKA material
index number does not overlap with one produced by FluDAG.**

The FluDAG calculation is now ready to run. First make a symbolic link from the
geometry file to a fixed file called ``dagmc.h5m``.
::

    $ ln -s geom.h5m dagmc.h5m

The problem can then be run with
::

    $ rfluka -e <path/to/fludag/executable/mainfludag> \
          ++{standard fluka options}++ <fludag_input_file>
