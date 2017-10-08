# dynaflow_output_processing

These folders contain matlab code for processing of the output files produced by [Dynaflow](https://blogs.princeton.edu/prevost/dynaflow/). They all handle output from a model of a hyperelastic material filled with trusses (dispersed state within the hyperelastic material) that are connected to the continuous hyperelastic phase through linear springs.

## TAPE88_process

This folder contains code to process the TAPE88 output files from Dynaflow. The y-axis can represent stress of the total material or other relationships of distances or forces between separate nodes. The x-axis are simulation time or strain of the total material.

The TAPE88 file contains data for particular nodes from each single time step.

## TAPE87_process

This folder contains code to process the TAPE87 and TAPE89 output files from Dynaflow. 2-D representations of the material can be produced for discreet times steps and .gif files can be produced that show a movie-like representation of the elongation process. In addition, the conductivity of the material can be plotted for the discreet time steps available in the output files.

The TAPE87 file contains data for all nodes for discrete time steps. The TAPE89 file contains data for elements for the same discreet times steps.
