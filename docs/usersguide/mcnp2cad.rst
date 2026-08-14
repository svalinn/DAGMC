mcnp2cad
========

UW has a tool, mcnp2cad_, which can
translate a MCNP model into ACIS format for use in future DAGMC simulations. The
tool builds CAD solids from MCNP cell descriptions often turning infinite bodies
and planes into finite versions. This tool has been integrated into Coreform Cubit,
and can be executed under the `File -> Import...` menu option by selecting the `MCNP`
file type.

.. figure:: assets/coreform_mcnp_import_screenshot.png

Executing this command will generate the MCNP model in the Coreform Cubit session
with material and cell importance assignment metadata included for addition to the
resulting DAGMC file upon export.

..  _mcnp2cad: https://github.com/svalinn/mcnp2cad
