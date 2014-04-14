geographic-eccentricity
=======================

Computes the geographic eccentricity and various other measures of 2-dimensional networks under various communication models. See http://arxiv.org/abs/1403.3007

Usage: geo_ecc [repeat x] [L x] [error x] \[truncate x\] ([rand x] | [sinr x y] | [exp x]) [file x]

This program invokes "voronoi", which is Steven Fortune's C program for computing 2-dimensional Voronoi diagrams, available at his webpage http://ect.bell-labs.com/who/sjf/
The voronoi program is invoked with the command 'system("./voronoi -t <points.tmp >delaunay.tmp");' The two files points.tmp and delaunay.tmp are thus created.
 
This program uses rand()/RAND_MAX to simulate a uniform random number distribution in [0,1]. It also uses the Box-Muller-Marsaglia transform to generate a Gaussian distribution from it. This may lead to the Neave anomaly depending on the C implementation (not tested).
 
Compile this program with '-lm' to link in the math library.
