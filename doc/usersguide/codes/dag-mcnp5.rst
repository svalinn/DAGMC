Code-specific steps for DAG-MCNP5
=================================

There are three varieties of code-specific steps for DAG-MCNP5:

1.  Defining attributes of the geometry using Cubit/Trelis groups
2.  Defining DAGMC runtime parameters using the DAGMC input file
3.  Specifying additional parameters on the command line

Geometry metadata
~~~~~~~~~~~~~~~~~

In DAG-MCNP5, the geometry file can be used to define material and density
assignments, and boundary conditions.

Materials and densities
-----------------------

The generic workflow description includes details on :ref:`grouping-basics`, but
a specific naming convention is required for DAG-MCNP5. To define materials,
both the MCNP material ID and density must be provided in the group name. The
format for the group name is as follows: ``mat:[matid]/rho:[density]``.

``[matid]`` should be replaced by the material ID that will be specified in the
MCNP input file. ``[density]`` should replaced by either the atomic density or
the mass density. As in the MCNP cell cards, positive values represent atomic
densities in [atoms/barn-cm] and negative values represent mass densities in
[g/cc].

For example, consider a problem where material 7 and an atomic density of 0.0223
should be assigned to volumes 4 through 18. The following command should be used
to specify that information:
::

    CUBIT> group "mat:7/rho:0.0223" add vol 4 to 18

All volumes must belong to a group, if you want volumes to be filled with vacuum
then add them to a group called, "mat:Vacuum":
::

    CUBIT> group "mat:Vacuum" add vol 4 to 18

If you would like to assign a material to the implicit complement, a special
procedure is needed. Since the implicit complement doesn't exist before running
DAGMC, and DAGMC can only recognize groups that contain an entity, the material
for the implicit complement must be specified as if it were being specified for
the graveyard volume. For example, if you would like the implicit complement to be
modeled as material 9 with a density of 1 g/cc, and the graveyard volume is
volume 102, the following command should be used:
::

    CUBIT> group "mat:9_comp/rho:-1_comp" add vol 102

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

    mat:Graveyard

For example, consider a geometry containing 99 volumes, all of which fit inside
a cube of side length 99 cm centered at the origin. The following commands would
create a valid graveyard for this problem:
::

    CUBIT> create brick x 100             # This will be volume 100
    CUBIT> create brick x 105             # This will be volume 101
    CUBIT> subtract vol 100 from vol 101  # This will produce volume 102
    CUBIT> group "mat:Graveyard" add vol 102

When DAG-MCNP5 is run, the importance any graveyard volumes will be set to zero.

It is still recommended that you create a graveyard volume even if your problem
has reflecting boundary conditions on all sides, although it is not strictly
necessary.

**Surface boundary conditions**

Surface boundary conditions can be specified for a given surface by adding the
surface to a group. The group names for reflecting and white boundary conditions
are ``boundary:Reflecting`` and ``boundary:White``, respectively. Note that periodic
boundary conditions are not currently supported.

For example, suppose you want to specify that surfaces 10 and 11 should be
reflecting surfaces. This command would achieve that:
::

    CUBIT> group "boundary:Reflecting" add surf 10 11

..  include:: dag-mcnp5_specific.txt


..  toctree::
    :hidden:

    dag-mcnp5_deprecated
