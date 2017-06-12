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
