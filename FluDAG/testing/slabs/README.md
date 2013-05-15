Introduction
================
This is the testing directory for the slabs benchmark. The benchmark consists of 10 cubes of side 10cm stacked in the z direction, with x-y axes centred along 0.0 in both directions, such that the total model (Non blackhole) extent is some 0.0 to 100.0 cm in the z direction and -5.0 to 5.0 in both x and y. The blackhole forms an annulus around the vacuum region and is 0.5 cm thick such that overall model dimensions are (-5.5,-5.5,0.5) to (5.5,5.5,100.5).

    *->

    +---+---+---+ ~~~~ +---+---+
    | 1 | 2 | 3 |      | 9 | 10|
    +---+---+---+ ~~~~ +---+---+

   0.0 cm                    100.0 cm

Source 
===============
The particle source used in 14.1 MeV neutrons where each volume is compsed entirely of vacuum. The source is distributed uniformly in x and y and direction is defined so as to prouduce a parallel beam of neutrons with x-y coordinates in the range -5 to 5 cm and z ~ 0.0 cm.

Scoring
==============
There is a track length based estimator in every volume region, but we set the volume which is used to normalise the result to one. Such that the track length estimator behaves as only the sum of interaction lengths within the volume. Since the volumes are composed of vacuum, and the neutron tracks behave as rays bounded by each cell, the maximum contribution any one ray contributes to the result is the thickness of the cell, i.e. 10 cm. We therefore define a successful calculation to return the score 10.0 for every cell.

Validation 
=============
The fiducial results as determined by running this benchmark are summarised below for both FLUKA and FluDAG

| Volume Number | Fluka | FluDAG |
|:-------------:|:-------------:|:-------------:|
    1    |   9.999996  |   9.999996
    2    |   9.999996  |   9.999996
    3    |   9.999996  |   9.999996
    4    |   9.999996  |   9.999996
    5    |   9.999996  |   9.999996
    6    |   9.999996  |   9.999996
    7    |   9.999996  |   9.999996
    8    |   9.999996  |   9.999996
    9    |   9.999996  |   9.999996
    10   |   9.999996  |   9.999996

These results are correct as 05/15/2013




 