
HTLMReportFileStructure
=======================

This project was designed to quickly provide an html report 
for the exploration of data gathered from Fluidigm’s biomark platform.

The input file format should be the "heat map" out put form the
biomark platform.

This analysis is heavily dependant on linear algebra and the reduction 
of a 96 dimensional gene space into a 1 dimensional progression axis.

The 96 genes are then analyzed against this progression vector defined 
as the difference between our start and end populations.

Single cells exposed to oskm are mapped into this one dimensional space
and the expression values inform a behavior of each of our genes along the
progression axis.

From a technical side this project involves Python R Bash and HTML,
Here I will only be pushing my own code and the pipeline in it’s entirety
will not at this time be available.

To see a demonstration :  ./runit.sh
Sentinel _   : Marks files to be read.
Sentinel __ : Marks files containing endpoint cell types. (skin and stem cells)
