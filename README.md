geographic-eccentricity
=======================

Computes the geographic eccentricity and various other measures of 2-dimensional networks under various communication models. See http://arxiv.org/abs/1403.3007

 Usage: [repeat x]   repeat experiment x times
        [L x]        simulates 4x² nodes scattered in a 4x * x rectangle
        [error x]    adds a position error with Gaussian distribution and
          standard deviation x
       [truncate x] removes links that are apparently longer than x
        [rand x]     uses the random graph model with link probability x
        [sinr x y]   uses the SINR model with min range x and max range y
        [exp x]      uses the expontential link model with average range x
        [file x]     appends compound results to the end of file x

This program invokes "voronoi", which is Steven Fortune's C program for computing 2-dimensional Voronoi diagrams, available at his webpage http://ect.bell-labs.com/who/sjf/
The voronoi program is invoked with the command 'system("./voronoi -t <points.tmp >delaunay.tmp");' The two files points.tmp and delaunay.tmp are thus created.
 
This program uses rand()/RAND_MAX to simulate a uniform random number distribution in [0,1]. It also uses the Box-Muller-Marsaglia transform to generate a Gaussian distribution from it. This may lead to the Neave anomaly depending on the C implementation (not tested).
 
Compile this program with '-lm' to link in the math library.
