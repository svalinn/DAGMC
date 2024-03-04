..  _h5m-file-spec:

DAGMC .h5m File Specification
=============================

The `.h5m` file is a generic file format supported by :term:`MOAB`. These files
may contain meshes intended for any number of purposes with various element
types. Some basics regarding :term:`MOAB` constructs is required to navigate
following specification, namely :term:`EntitySet`'s' and :term:`Tag`'s. See the
:ref:`glossary` for more information on these items.

Geometric EntitySets¶
---------------------

For an `.h5m` file to be used with :term:`DAGMC`, :term:`EntitySet`'s that are
"tagged" specific information must be present. Only :term:`EntitySet`'s that are
tagged with the required information will be "seen" by the :term:`DAGMC`
interface. These tags are used to identify the geometric entities (volumes,
surfaces, curves, and vertices) as well as their relationships to each other.
Tags on geometric :term:`EntitySet`'s are also used to establish topological
relationships. For example, the `GEOM_SENSE_2` tag (described in the
`geom_tags`_ table) is used to relate surfaces to the volumes on either side
of the surface.The `GEOM_DIM` tag is used to indicate the dimensionality of the
geometric `EntitySet`. See the `geom_dim_table`_ for valid entries for this tag.

Metdata EntitySets
------------------

To apply a :term:`DAGMC` geometry in transport, certain properties need to be
associated with the geometry. Examples of these properties include:

  - Material assignments
  - Boundary conditions
  - Temperatures
  - Tallies

Please refer to :ref:`code-specific-steps` for any properties that may be specific
to the transport code you intened to use.

.. _geom_tags:

.. table:: **Geometric EntitySet Tag Descriptions**

+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| Tag Name              | Type             | Real Type  | Size | Tagged On   | Purpose                                                                                                      |
+=======================+==================+============+======+=============+==============================================================================================================+
| `GLOBAL_ID`           | `MB_TYPE_INT`    | `int`      | 1    | `EntitySet` | Value of an ID associated with a geometric `EntitySet`.                                                      |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| `GEOM_SENSE_2`        | `EntityHandle`   | `uint64_t` | 2    | `EntitySet` | Relates a surface to the two volumes on either side of the surface. An entry in the first position           |
|                       |                  |            |      |             | indicates tht the surafce has a sense that is forward with respect to                                        |
|                       |                  |            |      |             | the volume `EntityHandle` in that position. An entry in the second position                                  |
|                       |                  |            |      |             | indicates that the surface has a sense reversed with respect to the volume `EntityHandle` in that position.  |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| `GEOM_SENSE_N_ENTS`   | `EntityHandle`   | `uint64_t` | N    | `EntitySet` | Relates a curve to any topologically adjacent surface `EntitySet`s.                                          |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| `GEOM_SENSE_N_SENSES` | `MB_TYPE_INT`    | `int`      | N    | `EntitySet` | Curve sense data correllated with the `GEOM_SENSE_N_ENTS` information.                                       |
|                       |                  |            |      |             | Values are `1` for a forward senses and `-1` for reversed senses.                                            |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| `CATEGORY`            | `MB_TYPE_OPAQUE` | `char`     | 32   | `EntitySet` | The geometric category of an `EntitySet`. One of "Vertex", "Curve", "Surface", "Volume", or "Group"          |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| `GEOM_DIM`            | `MB_TYPE_INT`    | `int`      | 1    | `EntitySet` | The dimensionality of a geometric `EntitySet`. See table below for meaning of values.                        |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+
| `NAME`                | `MB_TYPE_OPAQUE` | `char`     | 32   | `EntitySet` | A name assigned to an `EntitySet`. Use to indicate material assignments,                                     |
|                       |                  |            |      |             | boundary conditions, temperatures, and the implicit complement on                                            |
|                       |                  |            |      |             | `EntitySet`'s with a `CATEGORY` tag whose value is "Group"                                                   |
+-----------------------+------------------+------------+------+-------------+--------------------------------------------------------------------------------------------------------------+


.. _geom_dim_table:

.. table:: Dimensionality Values of the `GEOM_DIM` Tag

+-----------------+----------------------+
| Geometry Object | Dimensionality [*]_ |
+=================+======================+
| Vertex          | 0                    |
+-----------------+----------------------+
| Curve           | 1                    |
+-----------------+----------------------+
| Surface         | 2                    |
+-----------------+----------------------+
| Volume          | 3                    |
+-----------------+----------------------+

.. [*] The value of the `GEOM_DIM` tag on the geometric `EntitySet`.




.. table:: Oriented Bounding Box Tree Tag Descriptions

+------------+------------------+------------+------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------+
| Tag Name   | Type             | True Type  | Size | Purpose                                                                                                                                                              | Tagged On   |
+============+==================+============+======+======================================================================================================================================================================+=============+
| `OBB_ROOT` | `EntityHandle`   | `uint64_t` | 1    | This tag resides on geometric `EntitySet`'s. Its value is the handle of the associated OBB tree root `EntitySet`.                                                    | `EntitySet` |
+------------+------------------+------------+------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------+
| `OBB_GSET` | `EntityHandle`   | `uint64_t` | 1    | This tag resides on OBB tree root `EntitySet`'s. Its value is the handle of the associated geometric `EntitySet.`                                                    | `EntitySet` |
+------------+------------------+------------+------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------+
| `OBB`      | `MB_TYPE_DOUBLE` | `double`   | 9    | This tag resides on `EntitySets` in an OBB tree. The value of this tag is nine doubles representing the oriented bounding box for this `EntitySet`-node in the tree. | `EntitySet` |
+------------+------------------+------------+------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------+


Topology¶
~~~~~~~~

Every mesh-based geometry contains entity sets that are either
volumes, surfaces, or curves. There are two types of relationships that can
relate entities to other entities. The first is called a parent-child
relationship. Volumes are parents to surfaces that make up that volume; surfaces
are parents to curves; and curves are parents to the geometric vertices.

The second type of relationship is the set relationship, which is different from
a parent-child relationship. Each surface and curve is an entity set. The
surface entity sets contain the triangles and their vertices for that surface.
The curve entity sets contain edges and their vertices. The volume entity sets,
however, are empty. While a volume is parent to surfaces (the parent-child
relationship), the volume does not contain any mesh entities.

Sense tags¶
~~~~~~~~~~

Each surface is tagged with the two volume handles of the adjacent
volumes. The first of the two surfaces is designated as the forward direction
and the second is designated with the reverse direction. It is important to note
that these surfaces senses may not be consistent with how an MC code determines
the surface sense.