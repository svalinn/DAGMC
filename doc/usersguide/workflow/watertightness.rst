Watertightness
====================
The faceting from Cubit/Trelis is not a priori garenteed to be watertight, by which we 
define water tight to mean that the faceting (also known as tessellation) of two topologically
linked surfaces does not nessessarily have to be coincident. An example of the use of *make_watertight*
is show below, the black lines show the geometry after *make_watertight* the red shows the geometry before.

.. image:: watertight.png
   :height: 500
   :width:  600
   :alt: An example of the use of *make_wateright*, the black lines show the geometry after *make_watertight*
         the red shows the geometry before


make_watertight
---------------
A tool provided within the DAGMC distribution is the *make_wateright* tool. This tool takes the faceted
curve information and uses it to seal the triangle facets that meet on the same curve, thus enforcing
water tightness.

check_watertight
---------------
A tool provided within the DAGMC distribution is the *check_watertight* program, which will let the user
know of how watertight a given model is. It should be used as a post process to *make_watertight*, to 
verify that the model has been improved.
