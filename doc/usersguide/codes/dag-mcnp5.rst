Code-specific steps for DAG-MCNP5
=================================

There are three varieties of code-specific steps for DAG-MCNP5:

1. Defining attributes of the geometry using Cubit groups
2. Defining DAGMC runtime parameters using the DAGMC input file
3. Specifyinng additional parameters on the command line

Geometry metadata
~~~~~~~~~~~~~~~~~

In DAG-MCNP5, the geometry file can be used to define material and density
assignments, boundary conditions, and tallies.

Materials and densities
-----------------------

The generic workflow description includes details on :ref:`grouping-basics`, but
a specific naming convention is required for DAG-MCNP5. To define materials,
both the MCNP material ID and density must be provided in the group name. The
format for the group name is as follows: ``mat_[matid]_rho_[density]``.

``[matid]`` should be replaced by the material ID that will be specified in the
MCNP input file. ``[density]`` should replaced by either the atomic density or
the mass density. As in the MCNP cell cards, positive values represent atomic
densities in [atoms/barn-cm] and negative values represent mass densities in
[g/cc].

For example, consider a problem where material 7 and an atomic density of 0.0223
should be assigned to volumes 4 through 18. The following command should be used
to specify that information:
::

    CUBIT> group "mat_7_rho_0.0223" add vol 4 to 18

If a volume is not assigned to a material group, it will be treated as void by
DAGMC. This can be a fairly useful debugging tool when trying to identify
volumes that were not assigned to their appropriate group.

If you would like to assign a material to the implicit complement, a special
procedure is needed. Since the implicit complement doesn't exist before running
DAGMC, and DAGMC can only recognize groups that contain an entity, the material
for the implicit complement must be specified as if it were being specified for
the graveyard volume. For example, if you like the implicit complement to be
modeled as material 9 with a density of 1 g/cc, and the graveyard volume is
volume 102, the following command should be used:
::

    CUBIT> group "mat_9_rho_-1_comp" add vol 102

DAGMC will recognize that volume 102 is the graveyard, and the ``_comp`` keyword
will trigger it to assign the specified material and density to the implicit
complement rather than the graveyard.

Boundary conditions
-------------------

There are two general classes of boundary condition supported by DAG-MCNP5:
vacuum boundaries and surface boundary conditions.

**Vacuum boundaries: defining the "graveyard"**

Typical MCNP5 models contain a zero-importance volume that bounds the volumes of
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

To indicate to MCNP that a given volume is the graveyard volume, you must assign
it to a group with one of these names:
::

    graveyard
    outside.world
    rest.of.world

For example, consider a geometry containing 99 volumes, all of which fit inside
a cube of side length 99 cm centered at the origin. The following commands would
create a valid graveyard for this problem:
::

    CUBIT> create brick x 100             # This will be volume 100
    CUBIT> create brick x 105             # This will be volume 101
    CUBIT> subtract vol 100 from vol 101  # This will produce volume 102
    CUBIT> group "graveyard" add vol 102

When DAG-MCNP5 is run, the importance any graveyard volumes will be set to zero.

It is still recommended that you create a graveyard volume even if your problem
has reflecting boundary conditions on all sides, although it is not strictly
necessary.

**Surface boundary conditions**

Surface boundary conditions can be specified for a given surface by adding the
surface to a group. The group names for reflecting and white boundary conditions
are ``spec.reflect`` and ``white.reflect``, respectively. Note that periodic
boundary conditions are not currently supported.

For example, suppose you want to specify that surfaces 10 and 11 should be
reflecting surfaces. This command would achieve that:
::

    CUBIT> group "spec.reflect" add surf 10 11

Tallies
-------

You can also specify tallies in the geometry, although it is not required.
Tallies can be added by creating a group of volumes or surfaces and encoding the
metadata in the names of those groups. The group names should follow this
format:
::

    tally_[Cubit_tally_ID].[tally_type_keyword].[particles]

The ``[Cubit_tally_ID]`` field is an integer from 0 to 99.  Different tally
types may have the same Cubit ID and will still be consistent. The tally number
in MCNP will be 10 times the Cubit ID plus the tally type index (e.g. 4 for cell
flux tallies).

The ``[tally_type_keyword]`` should be one of the following:

+----------+------------------+
|Tally Type|tally type keyword|
+----------+------------------+
|f1        |surf.current      |
+----------+------------------+
|f2        |surf.flux         |
+----------+------------------+
|f4        |cell.flux         |
+----------+------------------+
|f6        |cell.heating      |
+----------+------------------+
|f7        |cell.fission      |
+----------+------------------+
|f8        |pulse.height      |
+----------+------------------+
|+f8       |qpulse.height     |
+----------+------------------+

It is possible to obtain tally results multiplied by particle energy (e.g. a
\*f2 tally in MCNP) by placing an ``e`` before the tally type. For example,
to make a ``*f2`` tally, the keyword should be ``esurf.flux``.

The ``[particles]`` tag is a string stating which particles will be
tallied.  To tally both photons and neutrons, set the tag to "np".
The default is neutrons only.  Should this be tag be omitted, only
neutrons will be tallied.

Here are some example Cubit commands to create tallies:
::

    CUBIT> group "tally_0.surf.current" add surf 1 to 4
    CUBIT> group "tally_0.cell.flux.p" add vol 7
    CUBIT> group "tally_1.ecell.heating.np" add vol 2 6
    CUBIT> group "tally_6.cell.heating.n" add vol 2 6
    CUBIT> group "tally_7.cell.flux.p" add vol 1 to 3
    CUBIT> group "tally_12.pulse.height.p" add vol 10 to 14
    CUBIT> group "tally_14.qpulse.height.p" add vol 10 to 14

The above are equivalent to following MCNP definitions:
::

    f1:n 1 2 3 4 T
    f4:p 7 T
    *f16:n,p 2 6 T
    f66:n 2 6 T
    f74:p 1 2 3 T
    f128:p 10 11 12 13 14 T
    +f148:p 10 11 12 13 14 T

Note that a total tally bin is always added.

.. _additional_parameters:

The DAG-MCNP5 input file
~~~~~~~~~~~~~~~~~~~~~~~~

DAG-MCNP5 input files should only contain the data card block of a standard
MCNP5 input file. There should be no cell or surface cards.

If you assigned materials to volumes using Cubit groups, you will need to define
these materials in the data cards. You will not need to do this if you are using
a UWUW material library.

If you did not assign tallies using Cubit groups, you will need to define them
in the data cards. If you did assign them using groups, you do not need to
define them in the data cards, although you can include tally modifiers here if
needed.

A new data card ``dagmc`` has been added to DAG-MCNP5 to define parameters for
the DAGMC geometry capability.
::

    Form: dagmc  keyword1=value   keyword2=value
           check_src_cell: behavior of CEL variable in SDEF card
                           on  [default] standard interpretation for
                                         CEL variable: source rejection
                           off           no cell rejection - assume that
                                         sampled position is in cell CEL
        overlap_thickness: allows particle tracking through small overlaps
                           {real} [default=0.0]
                   usecad: toggle usage of solid model geometry
                           off [default] ray-tracing limited to facets
                           on            ray-tracing performed on solid model
                                         geometry surfaces
                distlimit: toggle usage of flight distance sampled from
                           physics to accelerate ray-tracing search
                           off [default] do not use physics flight distance
                           on            do use physics flight distance

Running DAG-MCNP5
~~~~~~~~~~~~~~~~~

Running DAG-MCNP5 is identical to running the standard MCNP5, but a few new
keywords have been added to the command line to specify the necessary files.

:``gcad=<geom_file>`` (required): Specify the filename of the input geometry
    file. It can be in one of two formats: the MOAB (\*.h5m) format (this is the
    format produced by ``export_dagmc`` in Trelis/Cubit), or a facet file
    produced by DAGMC. If this entry is not present, DAG-MCNP5 will assume that
    it is running in standard MCNP5 mode.

:``fcad=<facet_file>`` (optional): [default: fcad] Specify the filename of the
    output facet file. This is the file produced by DAGMC that contains the
    geometry as well as the products of a number of preprocessing steps, which
    can be quite time-consuming. This file can be used as input with the
    ``gcad=`` keyword in subsequent runs to avoid spending time redoing the
    preprocessing steps.

:``lcad=<log_file>`` (optional): [default: lcad] Specify the filename of the
    output log file. This is a skeleton of an MCNP file which contains
    information about the volumes, surfaces, materials and material assignments,
    boundary conditions, and tallies defined in the geometry. If you specify a
    name other than the default for this file on the command-line, that file
    will be used instead of the one generated automatically by DAG-MCNP5. This\
    can be useful if you want to make small changes to your material
    assignments, importances, etc., but it cannot be used to change anything
    about the geometry. It is up to the user to ensure that the geometry file
    and log file being used correspond to each other. This runtime parameter is
    unique to DAG-MCNP5.
