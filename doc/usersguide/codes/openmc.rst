..  _OpenMC: https://mit-crpg.github.io/openmc

Code-Specific steps for OpenMC
==================================

**Note: DAGMC simulations in OpenMC also support the `UWUW <../uw2.html>`_
 workflow.**

There are two varieties of code-specific steps for OpenMC_:

1.  Defining attributes of the geometry using Cubit/Trelis groups
2.  Defining DAGMC runtime parameters using the OpenMC input files

Geometry metadata
~~~~~~~~~~~~~~~~~

The DAGMC geometry file can be used to define material assignments, boundary
conditions, and temperatures.

Materials and densities
-----------------------

The generic workflow description includes details on :ref:`grouping-basics`, but
a specific naming convention is required for DAG-OpenMC. To define materials,
the material ID must be provided in the group name. The format for the group
name is as follows: ``mat:[matid]``.

``[matid]`` should be replaced by the material ID that will be specified in the
``materials.xml`` file.

For example, consider a problem where material 7 should be assigned to volumes 4
through 18. The following command should be used to specify that information:
::

    CUBIT> group "mat:7" add vol 4 to 18

All volumes must belong to a group, if you want volumes to be filled with vacuum
then add them to a group called, "mat:Vacuum":
::

    CUBIT> group "mat:Vacuum" add vol 4 to 18

Boundary conditions
-------------------

**Surface boundary conditions**

The following boundary conditions available in OpenMC are supported in
DAGMC. These are:

- vacuum
- reflective
- transmission

For example, suppose you want to specify that surfaces 10 and 11 should be
reflecting surfaces. This command would achieve that:
::

    CUBIT> group "boundary:Reflecting" add surf 10 11

Surfaces without a specified boundary condition will be set to ``transmission``.

**Vacuum boundaries: defining the problem boundary**

The DAGMC model should have a "containing volume" which bounds the volumes of
interest. This volume should have surfaces that either remove particles from the
problem or reflect them back toward regions of interest. This volume should
surround the entire geometry with a shell of finite thickness. Any
geometric shape can be used, but a cubic shell is recommended. This
volume represents the boundary between the problem and the outside world.

To create a containing volume, create two volumes in Cubit with the same shape
and same center with one slightly larger than the other, making sure that both
bound the entire problem geometry. Then, subtract the smaller one from the
larger one. The result is a containing volume for the problem.

As mentioned above, this volume's surfaces should be assigned either vacuum or
reflective boundary conditions as is appropriate for the problem.

For example, consider a geometry containing 99 volumes, all of which fit inside
a cube of side length 99 cm centered at the origin. The following commands would
create a valid graveyard for this problem:
::

    CUBIT> create brick x 100             # This will be volume 100
    CUBIT> create brick x 105             # This will be volume 101
    CUBIT> subtract vol 100 from vol 101  # This will produce volume 102
    CUBIT> list volume 102 geom
    # collect surfaces of 102 from output #
    CUBIT> group "boundary:Vacuum" add <volume 102 surfaces>

Temperatures
~~~~~~~~~~~~

Cell temperatures can be defined in OpenMC using a similar syntax to boundary
conditions but with the "temp" keyword:
::
    CUBIT> group "temp:100" add vol x

All temperatures are assumed to be in Kelvin.

Implicit complement materials
-----------------------------

If you would like to assign a material to the implicit complement, a special
procedure is needed. Since the implicit complement doesn't exist before running
DAGMC, and DAGMC can only recognize groups that contain an entity, the material
for the implicit complement must be specified as if it were being specified for
the containing volume. For example, if you would like the implicit complement to
be modeled as material 9, and the containing volume is
volume 102, the following command should be used:
::

    CUBIT> group "mat:9_comp" add vol 102

DAGMC will recognize that volume 102 is the containing volume, and the ``_comp``
keyword will trigger it to assign the specified material and density to the
implicit complement rather than the containing volume.

Running DAG-OpenMC
~~~~~~~~~~~~~~~~~~

The command for running OpenMC is identical to an OpenMC run using native
geometry. Certain modifications to the OpenMC input files are required,
however. The element ``<dagmc>true</dagmc>`` must be present in the
``settings.xml`` file, and the DAGMC geometry must be named of symbolically
linked as ``dagmc.h5m``.

..  toctree::
    :hidden:

    dag-mcnp_deprecated
