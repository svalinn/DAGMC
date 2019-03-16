..  _OpenMC: https://openmc.readthedocs.io/en/latest/index.html

**Note: OpenMC supports the recommended** :ref:`UWUW`

|

Code-Specific steps for OpenMC
==================================

There are two varieties of code-specific steps for OpenMC_:

1.  Defining attributes of the geometry using Cubit/Trelis groups
2.  Defining DAGMC runtime parameters using the OpenMC input files

Geometry metadata
~~~~~~~~~~~~~~~~~

The DAGMC geometry file can be used to define material assignments, boundary
conditions, and temperatures.

Materials
---------

The generic workflow description includes details on :ref:`grouping-basics`, but
a specific naming convention is required for OpenMC. To define materials,
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


Assigning Materials by Name (recommended)
-----------------------------------------

OpenMC materials can also be assigned by name:
::

    CUBIT> group "mat:fuel" add vol 4 to 18

For this example, a material with the name attribute "fuel" must be present in
OpenMC's ``materials.xml`` file:
::
   <materials>
     <material id="40" name="fuel">
       <density units="g/cc" value="11" />
       <nuclide ao="1.0" name="U235" />
     </material>
   </materials>

This method for assigning materials is recommended for use with OpenMC as it
provides a more verbose description of the material definition and purpose.

This method of assignment also allows an easy transition to the UWUW workflow
for future models. One can embed a PyNE material library in the DAGMC model at
any point using ``uwuw_preproc`` without modification to the material
assignments to obtain a working UWUW model, provided that materials with
corresponding names and appropriate definitions are in the PyNE material
library.

**Note: material names must be unique in the materials.xml file for this style**
**of material assignment to work properly.**

Surface boundary conditions
----------------------------

The following boundary conditions available in OpenMC are supported in
DAGMC. These are:

- vacuum
- reflective
- transmission

For example, suppose you want to specify that surfaces 10 and 11 should be
reflecting surfaces. This command would achieve that:
::

    CUBIT> group "boundary:Reflecting" add surf 10 11

**Note: surfaces without a specified boundary condition will be set to** ``transmission`` **.**

Problem boundary
----------------

The DAGMC model should have a "containing volume" which bounds the volumes of
interest. This volume represents a particle "graveyard" or a region where
particles are removed from the simulation upon entry. This volume represents the
boundary between the problem and the outside world.  This volume should surround
the entire geometry with a shell of finite thickness. Any geometric shape can be
used, but a cubic shell is recommended to maximize performance.

To create a containing volume, make two volumes in Cubit/Trelis with the same
shape and same center with one slightly larger than the other, making sure that
both bound the entire problem geometry. Then, subtract the smaller one from the
larger one. The result is a containing volume for the problem.

For example, consider a geometry containing 99 volumes, all of which fit inside
a cube of side length 99 cm centered at the origin. The following commands would
create a valid graveyard for this problem:
::

    CUBIT> create brick x 100             # This will be volume 100
    CUBIT> create brick x 105             # This will be volume 101
    CUBIT> subtract vol 100 from vol 101  # This will produce volume 102
    CUBIT> group "mat:Graveyard" add volume 102

A volume in the ``mat:Graveyard`` group will be assigned a void material and its
surfaces will be given vacuum boundary conditions when the OpenMC simulation is
initialized.

Temperatures
------------

Volume temperatures can be defined in OpenMC using a similar syntax to materials or boundary
conditions but with "temp" as the keyword.

To assign a temperature of 900K to a volume one can use the following command.
::

    CUBIT> group "temp:900" add vol x

**Note: all temperatures are assumed to be in units of Kelvin in OpenMC.**

Implicit complement materials
-----------------------------

If you would like to assign a material to the implicit complement, a special
procedure is needed. Since the implicit complement doesn't exist before running
DAGMC and DAGMC can only recognize groups that contain an entity, the material
for the implicit complement must be specified as if it were being specified for
the graveyard volume. For example, if you would like the implicit complement to
be modeled as material 9, and the graveyard volume is volume 102, the following
command should be used:
::

    CUBIT> group "mat:9_comp" add vol 102

DAGMC will recognize that volume 102 is the graveyard volume, and the ``_comp``
keyword will trigger it to assign the specified material and density to the
implicit complement rather than the containing volume.

Running DAG-OpenMC
~~~~~~~~~~~~~~~~~~

The command for running OpenMC is identical to an OpenMC run using native
geometry. Certain modifications to the OpenMC input files are required,
however. The element ``<dagmc>true</dagmc>`` must be present in the
``settings.xml`` file, and the DAGMC geometry must be named or symbolically
linked as ``dagmc.h5m`` in the directory where the ``openmc`` command is
executed.

..  toctree::
    :hidden:

    dag-mcnp_deprecated
