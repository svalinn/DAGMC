mcnp2cad
========

UW has a tool, mcnp2cad_, which can
translate a MCNP model into ACIS format for use in future DAGMC simulations. The
tool builds CAD solids from MCNP cell descriptions often turning infinite bodies
and planes into finite versions. At the time of writing, the only unsupported
MCNP surface descriptions are limited to GQ's and SQ's. To run mcnp2cad all that
is needed is an MCNP input deck,
::

    $ mcnp2cad test.inp

Will result in a file called out.sat which will contain the CAD version of your
MCNP input. Furthermore, mcnp2cad automatically transfers material and
importance assignments into the CAD model and will be translated to the DAGMC
file when processed.

..  _mcnp2cad: https://github.com/svalinn/mcnp2cad
