# Gibbs-Sampling-for-Multiple-Sequence-Alignment
Written in C after Lawrence et al. 1993 in Science. I wrote this algorithm for a masters course in bioinformatics at Sorbonne Université in Paris. When given a Fasta file with a list of sequences (from any alphabet--DNA, RNA, amino acid, whatever you'd like) it will find subtle similarities between those sequences, even if the similarities are not exact.
The algorithm runs in O(m*n), where m is the number of sequences inputted and n is the length of each sequence. 

The following is a description of each file in the repository:

Fasta files for testing the algorithm:

  1. TestSeq2.fasta
    A challenging edge-case for the algorithm, which tends to get stuck on the ends of these sequences unless it can phase shift. 
  2. TestSeq3.fasta
    A basic set of three sequences for testing.
  3. lipocalin.fst
    A set of 5 actual protein sequences. The results of the algorithm run on this file for motif lengths of 6 and 10 can be found in Gibbs Sampling Results.pdf.
  
Files containing the functions and main for the algorithm:

  1. gibbs_main.c
    The main function for the algorithm
  2. gibbs_functions.h
    All helper functions that need to be called in the main function. 
  
Results Files:

  1. Gibbs Sampling Results.pdf
    Simple output from the most recent version of the algorithm, showing results from the lipocalin.fst Fasta file for motif lengths of 6 and 10. 
  2. Papier sur Gibbs Sampling.pdf
    A more in-depth explanation of the algorith (albeit on an older, buggier version). This file is unfortunately in French, but if it can be read it is quite informative. 
    
For all additional information on the algorithm, see Lawrence et al. 1993 (10.1126/science.8211139). I actually studied with Charles Lawrence at the Center for Computational Molecular Biology at Brown, so it was a pleasure to re-construct his most famous publication. 
  
