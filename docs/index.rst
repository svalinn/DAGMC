.. Double-Down documentation master file, created by
   sphinx-quickstart on Fri Mar 26 13:29:08 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DAGMC:  The Direct Accelerated Geometry Monte Carlo Toolkit
===========================================================

DAGMC is a toolkit for performing Monte Carlo radiation transport simulations on
CAD-based geometry models. It is built on top of the MOAB mesh database. DAGMC
is developed and maintained by the Computationan Nuclear Engineering Research
Group (CNERG) at the University of Wisconsin -- Madison.

.. carousel::
    :show_controls:
    :show_fade:
    :show_captions_below:

    .. figure:: assets/manifold-cad.png
       :height: 400px

       CAD

       A pipe manifold modeled in :term:`Coreform Cubit`.

    .. figure:: assets/manifold-tris.png
       :height: 400px

       Surface Mesh w/ Embedded Topology

       A DAGMC surface mesh of the piping manifold.

    .. figure:: assets/manifold-flux.png
       :height: 400px

       Flux Results

       A flux mapping generated using :term:`OpenMC`.

Usage of CAD-Based Geometry vs. CSG in Monte Carlo Particle Transport
=====================================================================

Most Monte Carlo codes natively support :term:`CSG` geometry as it is very robust
for particle tracking applications. However, use of  :term:`CAD` (Computer-Aided
Design) geometry offers several advantages, especially in the context of modern
engineering and design practices. :term:`CAD` geometry provides benefits in a
number pieces in the engineering design chain:

Visualization and Realism
-------------------------
:term:`CAD` geometry provides an interacctive visual representation of the
design, making it easier for designers and engineers to conceptualize and
iterate upon their ideas. Unlike :term:`CSG`, which relies on a fixed set of
surface objects, :term:`CAD` geometry offers a more intuitive approach by
allowing users to see the design in a realistic manner.

Flexibility and Iteration
--------------------------
:term:`CAD` geometry offers unparalleled flexibility in design iteration.
Designers can easily modify shapes, dimensions, and features with simple clicks
and adjustments, enabling rapid prototyping and experimentation. In contrast,
while most Monte Carlo codes :term:`CSG` offer some option to visualize geometry,
these tools are often limited (2D slices and laggy interaction) and require a
specific native format. Iteration on :term:`CSG` geometry is in turn inherently
less intuitive.

Simulation and Analysis
------------------------
:term:`CAD` geometry enables integration with simulation and analysis tools for
evaluating the performance and behavior of designs under various conditions.
From :term:`FEA` to fluid dynamics simulations, :term:`CAD` software allows
engineers to validate their designs before physical prototyping, saving time and
resources.

One significant advantage of CAD-based geometry in simulation and analysis is its
capability for multiphysics domain mapping. Engineers can model complex systems
involving multiple physical phenomena, such as structural mechanics, heat
transfer, and fluid flow, by seamlessly integrating different simulation modules
within a CAD environment. This allows for a comprehensive understanding of how
different aspects of the design interact with each other, leading to more
accurate predictions and optimized designs.

CAD geometry also facilitates the mapping of simulation results back to the
design, providing valuable insights for further refinement. Engineers can
identify critical areas of stress, temperature gradients, or fluid flow
restrictions directly on the CAD model, enabling targeted design improvements.

While CSG can theoretically support similar analyses, CAD geometry offers a more
practical and integrated approach, simplifying the workflow and enhancing
productivity.

Advantages and Disadvantages of Surface Mesh (Triangles) vs. Volumetric Mesh (Tetrahedra)
===========================================================================================

Advantages of Surface Mesh (Triangles)
----------------------------------------

2. **Higher Fidelty Boundary Representations**: Volumetric meshes are often
   limited in how well they can resolve the boundary between parts in a
   :term:`CAD` model due to constriants on mesh quality of interior elements.
   This commonnly results in a more coarse approximation of boundaries than can
   be achieved with a surface mesh for the same number of triangle elements.
   They are also able to more accurately capture features of varying sizes, such
   as sharp corners, thin features, or regions of high curvature. This in turn
   translates to more accurate representation of surface area and volume.

1. **Robust Meshing**: A surface mesh can be generated for nearly any manifold
   volume. It is very rare that a triangle mesh cannot be generated for a given
   volume.

1. **Fewer Boundary Crossings**: For the same level of detail, surface meshes
   typically have fewer boundary crossings compared to volumetric meshes. This
   can reduce the computational overhead associated with tracking particles
   through volumetric mesh elements, making surface meshes more efficient for
   certain regions in which the particle's average path length is much larger
   than the local elements.

Disadvantages of Surface Mesh (Triangles)
-------------------------------------------

1. **Limited Volume Representation**: Surface meshes do not capture the internal
   volume of the object directly, which can limit their applicability for
   simulations involving volumetric phenomena.

2. **Ray Tracing Operations**: While particles are able to travel from one
   surface to another in a volume containing a low- or zero-density material,
   the ray tracing operations involved in this process can be computationally
   expensive compared to adjacency search when traversing volumetric elements.

3. **Less Accurate Volume-based Results**: Since surface meshes do not directly
   model the volume of the object, simulations relying on volumetric quantities
   may yield less accurate results compared to volumetric meshes. Extrapolating
   volume-based information from surface meshes can introduce uncertainties and
   approximation errors.

Advantages of Volumetric Mesh (Tetrahedra)
------------------------------------------

1. **Volume Representation**: Volumetric meshes directly represent the internal
   volume of the object, allowing for a more accurate simulation of complex
   three-dimensional phenomena and property representation (density, temperature, etc.)
   during particle transport.

2. **Accurate Boundary Representation**: With a volumetric mesh, the surface
   geometry is implicitly defined by the tetrahedral elements, ensuring accurate
   representation of the boundaries and interfaces between different materials
   or regions.

Disadvantages of Volumetric Mesh (Tetrahedra)
-----------------------------------------------

1. **Higher Computational Cost**: Generating and solving volumetric meshes can
   be computationally expensive, especially for large and complex geometries.
   The presence of tetrahedral elements throughout the volume increases the
   number of degrees of freedom and requires more computational resources.

2. **Mesh Quality Concerns**: Ensuring high-quality tetrahedral meshes, such as
   avoiding element distortion or ensuring element aspect ratios, can be
   challenging, particularly for irregular or highly curved geometries. Poor
   mesh quality can cause particle tracking problems. Additionally, the set of
   models that can be represented by a volume mesh is more limited than the set of
   models that can be represented by a surface mesh.

In summary, the choice between volumetric mesh (tetrahedra) and surface mesh
(triangles) depends on the specific requirements of the simulation or analysis,
including the desired level of accuracy, computational resources, and the nature
of the geometry being modeled. Volumetric meshes offer accurate property
representation across a single volume/part and are suitable for tallying
internal fields while surface meshes are generally able to better conserve
volume and surface area and may be better suited for models involving complex
geometry or large numbers of parts. Performance considerations between the two
are often problem-dependent and should be evaluated on a case-by-case basis.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Table of Contents
=================

.. toctree::
    :maxdepth: 2

    self
    usersguide/index
    developerguide/index
    theoryguide/index
    filespec/index
    api/index
    glossary
