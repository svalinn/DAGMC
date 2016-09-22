Frequently Asked Questions
==========================

* Q1_ How does the implicit complement work?
* Q2_ Can I automate the CUBIT conversion process?
* Q3_ The install procedure is tedious, isn't there another way?
* Q4_ My problem doesnt work, it ...

.. _Q1:

Q1: How does the implicit complement work?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**A1** Since your geometry has been imprinted and merged, all surfaces fall into one of two categories:

1. surfaces that are shared by two volumes
2. surfaces that are used in only one volume

All surfaces in this second group are, by definition, on the boundary
of the complement region since they have an explicitly defined region
on one side and nothing on the other.  The implicit complement is
formed automatically from this list of surfaces.

From the point of view of MCNP, and extra cell is created to represent
the implicit complement.  Your output file will include one entry for
each volume in your geometry, including the graveyard, plus one extra
entry representing the complement.

.. _Q2:

Q2: The CUBIT conversion process is tedious and frustrating. Is there a way to avoid or automate most of this work?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**A2** Yes, a script has been written that will automate all the
necessary steps to perform the CUBIT conversion process with limited
user input. Information about this script can be found on the
AutomatedCubitConversion website.

.. _Q3:

Q3: The install procedure is tedious, isn't there another way?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**A3** Yes, we have written an install script that simplifies the building of DAGMC and its subsequent programs, we
recongnize that the dependency stack of DAGMC is large, and a pain to install, but the benefits are that very complex
geometric problems are now tractable. We also plan on adding more build options in the future.

.. _Q4:

Q4: My problem doesnt work, it crashes, loses particles, etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**A4** Occasionally users experience problems using DAGMC, new workflows, different models, there are several potential
reasons why your specific problem is not behaving as you expect, we recommend any user that downloads DAGMC also joins
our `Google group <https://groups.google.com/forum/#!forum/dagmc-users>`_, its an excellent resource and also a repository
of answers to users previous issues.
