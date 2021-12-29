# Forward Wright-Fisher Simulator
The program given below is based on the algorithm described in Padhukasahasram et al. 2008.



http://www.sourceforge.net/projects/forwsim/files/forwsim.c

http://www.sourceforge.net/projects/forwsim/files/mtrand.h

http://www.sourceforge.net/projects/forwsim/files/mtrand.cpp



Download the files in the 3 links above. forwsim.c is a program that can efficiently simulate the evolution of a diploid Wright-Fisher population of constant size with uniform crossing-over rate, uniform mutation rate and with or without self-fertilization, forward in time. To compile type: g++ -O3 -o forwsim forwsim.c mtrand.cpp -Wno-deprecated. There are 8 different command line parameters which are:

-samples (Number of samples to output from final population.)

-pop (Total number of chromsomes in the diploid population. Should be a even number.)

-len (Length of the DNA sequence. Choose larger than del x u x pop)

-r (Per-generation per-sequence rate of recombination)

-u (Per-generation per-sequence rate of mutation)

-s (Probability of self-fertilization)

-gen (Total number of generations)

-del (Generations after which fixed mutations are removed)



To run type ./forwsim along with all the command line parameters in the same order. For example: ./forwsim -samples 100 -pop 1000 -len 1000000 -r 0.50 -u 0.50 -s 0.0 -gen 10000 -del 500. The output file is called finalpopulation.txt. Each line in this file corresponds to a chromosome in the population and the numbers represent positions of mutations in it. Consecutive pairs of lines represent homologous chromosomes in the diploid population.
