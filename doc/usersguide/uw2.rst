The University of Wisconsin Unified Workflow
----------------------------------------

The University of Wisconsin Unified Workflow (UW<sup>2</sup>) aims to solve the 
issue of running the same Monte Carlo problem using mutiple physics codes. Currently,
if you wish to run the same problem in multiple codes you must fully recreate the
input deck for each code, or maybe write a full syntax translator. UW<sup>2</sup>, allows users
to tag or associate groups of volumes or surfaces with a simple human readable syntax
that is translated and stored in the geometry file of a DAGMC problem.

The workflow uses the Python for Nuclear Engineering toolkit `PyNE <http://pyne.io>`_, we 
levereage the existing infrastructure in PyNE to allow a consistent transport problem to be
defined across all MC codes.

Materials
+++++++++++++++++++++++++++++++++++++++

Materials are the most painful transfer from code to code, since each MC code 
specifies materials in a different way. Instead, we tag groups of volumes
with a name and syntax that corresponds to material compositions in a predefined
material library.

The syntax for describing materials is:
::
     %> group "mat:<Name of Material>"

If you wish to specify a material that is present in the library with a different 
density, then 
::
     %> group "mat:<Name of Material>/rho:<density>"

So for example, to specify a Stainless Steel at a density of 12.0 g/cc,
::
     %> group "mat:Stainless Steel/rho:12.0"


Scoring
+++++++

Each MC code implements tallies or scores, in very specific ways such that there
is sometimes no equivlent to a tally you may be familiar with, code to code. However
, there is a syntax to allow you to request scores on geomemtric elments, for example,
::
     %> group "tally:Flux/Neutron"

or,
::
     %> group "tally:Photon/Current"

Using the underlying PyNE libraries, we can write out the appropriate MC code tally specification
snippet, this allows the number of codes the DAGMC supports to grow organically with those that
PyNE supports. When PYNE cannot fulfill your tally request it will warn you.

UWUW Data
+++++
The UWUW data is incorporated into the geometry file (*.h5m) file using a Python script, uwuw_preproc, 
the purpose of which is to take the users material library e.g. my_nuc_library.h5 and extract the materials
requested, placing them into the geometry file. Having already marked up your geometry using the methods
mentioned in previous sections, we can run the preprocess script,
::
   %> uwuw_preproc -f <dagmc h5m filename> -d <path to nuclear data library> \
                   -o <output h5m filename>

Be sure to examine the output of this script which will inform you of the materials and densities requested and 
also the list of tallies that were produced. A sample output is shown below
::
   %> uwuw_preproc -f test_geom.h5m -d $HOME/.local/lib/python2.7/site-packages\
                     /pyne/nuc_data.h5 -o output.h5m

Also, the script will fatal error if the material is not found in the material library
::
   %>uwuw_preproc -f test_geom.h5m -d $HOME/.local/lib/python2.7/site-packages \
                     /pyne/nuc_data.h5 -o output.h5m
   ...
   mat:Lead
   mat:Lead/rho:12.8
   mat:Beryllium
   mat:Tungsten
   mat:Graveyard
   mat:StainlessSteel
   Material {StainlessSteel} doesn't exist in pyne material lib

Gotchas
=======
When using the "-o" option with a filename that does not match the "-f" filename will produce a new file with all the material and tally data, but 
absent of the original geometry data. This feature allows the user to ensure that the uwuw_preproc script runs to succesful completion before 
being used on the original file. Once you have ensured a sucessful run, it is recommended that you run once more with the "-o" option set to the 
original filename i.e.
::
   %>uwuw_preproc -f test_geom.h5m -d $HOME/.local/lib/python2.7/site-packages \
                    /pyne/nuc_data.h5 -o output.h5m
   Success!!
   %>uwuw_preproc -f test_geom.h5m -d $HOME/.local/lib/python2.7/site-packages \
                    /pyne/nuc_data.h5 -o test_geom.h5m
  
The reason for this behaviour is because it can take some time to produce a workflow ready facet file, having done dagmc_preproc and then make_watertight
and so on.

Worked Example
+++++

Open Cubit, and lets place some volumes, create our first cube, we will create 4 cubes of side 10 cm, shifting each in a different direction
::
   %>brick x 10
   %>move Volume 1 x 20 include_merged
   %>group "mat:Lead" add volume 1
   %>group "tally:Photon/Flux" add volume 1
   %>brick x 10
   %>move Volume 2 x -20 include_merged
   %>group "mat:Lead" add volume 2
   %>group "tally:Photon/Flux" add volume 2
   %>brick x 10
   %>move Volume 3 y -20 include_merged
   %>group "mat:Lead/rho:12.3" add volume 3
   %>group "tally:Photon/Flux" add volume 3
   %>brick x 10
   %>move Volume 4 y 20 include_merged
   %>group "mat:Lead/rho:12.3" add volume 4
   %>group "tally:Photon/Flux" add volume 4
   %>brick x 100
   %>brick x 105
   %>subtract volume 5 from volume 6
   %>group "mat:Graveyard" add volume 7
   %>imprint body all
   %>merge all
   %>set attribute on
   %>export acis "example.sat" overwrite

Now the file is ready for preprocessing, first we must facet the file;
::
   %>dagmc_preproc example.sat -o example.h5m

Now we can insert all the material data we need;
::
   %>uwuw_preproc -f example.h5m -d $HOME/.local/lib/python2.7/site-packages\
                     /pyne/nuc_data.h5 -o example.h5m

Your output from this step should look exactly the same as below
::
   Making nuc_data at example.h5m
   skipping atomic mass data table creation; already exists.
   Materials Requested....
   mat:Graveyard
   mat:Lead
   mat:Lead/rho:12.3
   Tallies Requested....
   Photon Flux PHFLUX1
   Photon Flux PHFLUX2
   Photon Flux PHFLUX3
   Photon Flux PHFLUX4

So we see echoed back to us that we requested a Graveyard, and two different material assignments, one for Lead, 
as defined in the material library and another kind of Lead at a different density than the library version. We 
also see that 4 tallies were requested, the photon flux in each volume.

Example Input
======
We are now ready to run once we have made the input deck for each Monte Carlo code, we wish to launch 10^5 particles, 
from a point source located at 0 0 0, with isotropic angular behaviour with photons of 1 MeV. The input for MCNP and
FLUKA are shown below, MCNP for example let us call this mcnp.inp ;
::
   example of UWUW
   c notice no cell cards
   c notice no surface cards
   c notice no blank lines!
   sdef x=0.0 y=0.0 z=0.0 par=2 erg=1.0
   c notice no materials
   c notice no tallies
   mode p
   nps 1e5
   print 

And Fluka, let us called this fluka.inp;
::
   TITLE
   * Set the defaults for precision simulations
   DEFAULTS                                                              PRECISIO
   * Define the beam characteristics
   BEAM          -0.001             10000.0                              PHOTON
   * Define the beam position
   BEAMPOS           0.        0.        0.
   * Notice the FLUGG section
   GEOBEGIN                                                              FLUGG
   GEOEND
   * notice no material assignments
   * notice no scoring assignments
   * ..+....1....+....2....+....3....+....4....+....5....+....6....+....7...
   RANDOMIZ         1.0
   * Set the number of primary histories to be simulated in the run
   EMF
   START           1.E5
   STOP

MCNP Run
========
So we are now ready to run the example, first DAG-MCNP5;
::
   %> mcnp5 i=mcnp.inp g=example.h5m

You should see the following on screen
::
   The implicit complement's total surface area = 128550
   This problem is using DAGMC version    1.000 w/ DagMC r   0
   Using default writer WriteHDF5 for file fcad 
   /mnt/data/prod/uwuw_example/web_example/example.h5m
   Materials present in the h5m file
   mat:Lead
   mat:Lead/rho:12.3
   Tallies present in the h5m file
   PHFLUX1
   PHFLUX2
   PHFLUX3
   PHFLUX4
   Going to write an lcad file = lcad
   Tallies
             Thread Name & Version = MCNP5, 1.60
             Copyright LANS/LANL/DOE - see output file
                                     _                                      
               ._ _    _  ._   ._   |_                                      
               | | |  (_  | |  |_)   _)                                     
                               |                                            
           
   comment.  photon   importances have been set equal to 1.                                                               
   comment. using random number generator  1, initial seed = 19073486328125      
   Turned OFF ray firing on full CAD model.
   Set overlap thickness = 0
   imcn   is done
  
    warning.  material        1 has been set to a conductor.                                                               
    warning.  material        2 has been set to a conductor.                                                               
  
                              ctm =        0.00   nrn =                 0
   dump    1 on file runtpe   nps =           0   coll =                0
     xact   is done

   cp0 =   0.01
   run terminated when      100000  particle histories were done.
  
                                ctm =        0.05   nrn =            900033
   dump    2 on file runtpe   nps =      100000   coll =            56221
   mcrun  is done

Feel free to examine the output of the run, but this provides a simple example on what to
expect.

FluDAG Run
==========
And now FluDAG, first we produce the mat.inp snippet file, this must then be pasted into
the full Fluka input deck
::
   %> mainfludag example.h5m

The mat.inp file should look like
::
   *...+....1....+....2....+....3....+....4....+....5....+....6....+....7...
   ASSIGNMA       LEAD1        1.
   ASSIGNMA       LEAD1        2.
   ASSIGNMA       LEAD2        3.
   ASSIGNMA       LEAD2        4.
   ASSIGNMA    BLCKHOLE        5.
   ASSIGNMA      VACUUM        6.
   *...+....1....+....2....+....3....+....4....+....5....+....6....+....7...
   MATERIAL         82.   207.217     11.35       26.                    LEAD1     
   MATERIAL         82.   207.217      12.3       27.                    LEAD2     
   *...+....1....+....2....+....3....+....4....+....5....+....6....+....7...
   * UW**2 tallies
   * PHFLUX1
   USRTRACK         1.0    PHOTON       -21        1.1.0000e+03     1000.PHFLUX1
   USRTRACK       10.E1     1.E-3                                               &
   * PHFLUX2
   USRTRACK         1.0    PHOTON       -21        2.1.0000e+03     1000.PHFLUX2
   USRTRACK       10.E1     1.E-3                                               &
   * PHFLUX3
   USRTRACK         1.0    PHOTON       -21        3.1.0000e+03     1000.PHFLUX3
   USRTRACK       10.E1     1.E-3                                               &
   * PHFLUX4
   USRTRACK         1.0    PHOTON       -21        4.1.0000e+03     1000.PHFLUX4
   USRTRACK       10.E1     1.E-3                                               &

As of the current time you will need to add two lines manually, this is because the 
component of the code which identifies neutron cross section data is not yet complete.
::
   *...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....
   LOW-MAT        LEAD1       82.       -2.      296.                    LEAD 
   LOW-MAT        LEAD2       82.       -2.      296.                    LEAD 

This must be pasted into the Fluka input and then run, as you would any Fluka, but
with the exception that we give the rfluka script an exectuable argument, and  new
"-d" argument which specifies the geometry filename
::
   %> $FLUPRO/flutil/rfluka -N0 -M1 -e mainfludag -d example.h5m fluka.inp

The code should run and sucessfully produce the following
::
   $TARGET_MACHINE = Linux
   $FLUPRO = /mnt/data/opt/fluka/fluka/

   Initial seed already existing
   Running fluka in /mnt/data/prod/uwuw_example/web_example/fluka_26362

   ======================= Running FLUKA for cycle # 1 =======================

   Removing links
   Removing temporary files
   Saving output and random number seed
   Saving additional files generated
   Moving fort.21 to /mnt/data/prod/uwuw_example/web_example/fluka001_fort.21
   End of FLUKA run

DagSolid Run
============
DagSolid is probably the most trivial of all the UWUW enabled codes to run, copy the vis.mac file from DAGMC/geant4/build/vis.mac
::
   %> DagGeant4 example.h5m

After some loading you should see a GUI window open (if you build geant4 with visualisation on), we can then use the Geant4 general particle
source to emulate the behaviour of the previous two codes,
::
   Idle> /gps/particle gamma
   Idle> /gps/ang/type iso
   Idle> /gps/energy 1.0 MeV

Now we are ready to run,
::
   Idle> /run/beamOn 1000000


