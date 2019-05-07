#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>


//D'abord, on crée des nouveaux types, qui s'appellent tgSeq et tgMotif.
//First, we create two new types called tgSeq and tgMotif. 

//tgSeq est une collection d'information sur une seule séquence, qu'on peut
	//utiliser pour créer un tableau de séquences
//tgSeq is a collection of information on a single sequence that was read into the algorithm.
	//A table of tgSeq's is a table of sequences. 
typedef struct _tgSeq { char *titre ; char *seq ; int lg ; } tgSeq ;

/*
Our sequences might look like this:
1. B	C	A	D	C	B	D	A
2. B	D	C	B	C	B	A	D	B										
3. A	D	B	C	B	D	D	A	

For each of the example tables in the function descriptions, we're pretending that the algorithm
	has decided the following:	
	
	Sequence 3 is the excluded sequence
	
	The positions of the motifs in sequences 1 and 2 are indices 0 and 2, respectively,
		so the current alignment is like this:

	1:     [B C A D] C B D A 
	2: B D [C B C B] A D B 
*/

//tgMotif est une collection d'information sur une seule motif, qu'on peut manipuler 
	//pour éventuellement trouver le cadre séquence. 
//tgMotif is a collection of information on a single motif, which we can manipulate to 
	//eventually find the target sequence. 
typedef struct _tgMotif { int lgMotif ; int nbSeq ; tgSeq *tSeq ; int *tDebut ; } tgMotif;

//Ici, on défine le nombre de résidues. C'est important car on peut utiliser cet algorithme 
	//pour les séquences de protéines et aussi de ADN.
//Here, we define the number of possible residues (also known as the length of the alphabet).
	//It's important that we can use this algorithm on amino-acid sequences as well as DNA
	//and RNA.
#define NUM_RESIDUES 26
#define NUM_ITERATIONS 10000
#define NUM_COMMANDS 2
#define MAX_SEQ_LN 1000
#define LG_MOTIF 10
//Si PHASE_SHIFT_FREQUENCY > NUM_ITERATIONS, l'algorithme ne fait pas de phase shifts. 
//If PHASE_SHIFT_FREQUENCE > NUM_ITERATIONS, the algorithm will not perform phase shifts. 
#define PHASE_SHIFT_FREQUENCY 100
#define MAX_SHIFT 5


/*-----------------------------------------------------------------------------------------------*/
/*											readFasta											 */
/*-----------------------------------------------------------------------------------------------*/	

tgSeq* readFasta(char *nomFi, int *pNbSeq)
{
//This function takes in a fasta file and returns a table of tgSeq types. In other words, 
	//it converts a fasta file into a table of sequences that the algorithm can manipulate. 

	//Inputs: 
		//1. a pointer to char nomFi, a character name for the fasta file being loaded in
		//2. a pointer to int pNbSeq, whose value will also be set by this function. 

	//Outputs:
		//1. a table of tgSeq types, which is a table of the sequences we're dealing with

	//Values changed:
		//1. int pNbSeq, which will set to the number of sequences read from the fasta file. 
			

    tgSeq *tS=NULL;

    char ligne[MAX_SEQ_LN];
    int lgBuff;
    int nbSeq=0;
    FILE *pFi;

    pFi=fopen(nomFi, "rt"); /*Ouvrir le fichier (open the file) */
    if(pFi==NULL){
        fprintf(stderr, "readFasta:: fichier introuvable: %s\n", nomFi);
        return NULL;
    }

    /*Lire les  lignes (Read the lines)*/
    while(fgets(ligne, MAX_SEQ_LN, pFi)!=NULL) {
    	/*Calculer la longueur de la ligne*/
    	/*Calculate the length of the line*/
        lgBuff=strlen(ligne); 
        /*Enlever le retour à la ligne s'il y en a un*/
        /*Remove the return at the end of the line if there is one*/
        if(ligne[lgBuff-1]=='\n'){ 
            ligne[lgBuff-1]='\0';
            lgBuff--;
        }
        /*Si la ligne n'est pas vide*/
        /*If the line is not empty*/
        if (lgBuff>0) { 
            /*Si c'est une ligne de titre*/
            /*If it is a title line for a sequence*/
            if(ligne[0]=='>') { 
                nbSeq++;               
                /*Ajouter une case au tableau des séquences*/
                /*Add an element to the table of sequence*/
                tS=(tgSeq*)realloc(tS, sizeof(tgSeq)*(nbSeq));               
                if(tS==NULL){
                    fprintf(stderr, "readFasta:: Erreur Allocation tS\n");
                    exit(1);
                }
                /*créer les membres de la nouvelle séquence*/
                /*Create the members of the new sequence*/
                tS[nbSeq-1].lg=0;
                tS[nbSeq-1].seq=NULL;
                /*Allouer une chaine de caracteres pour le titre */
                /*Allocate a chain of characters for the title*/
                tS[nbSeq-1].titre=(char*)malloc(sizeof(char)*(lgBuff+1)); 
                if(tS[nbSeq-1].titre==NULL){
                    fprintf(stderr, "readFasta:: Pb Allocation titre\n");
                    exit(1);
                }
                strcpy(tS[nbSeq-1].titre, ligne); /*Recopier la ligne de titre*/ 
            /*C'est une ligne de séquence*/
            /*If it is a line of characters for the sequence*/
            } else {
                if(nbSeq==0){
                    fprintf(stderr, "readFasta:: format incorrect, exiting\n");
                    exit(1);
                }
                /*Allouer ou agrandir la chaine de caracteres pour la sequence */
                /*Allocate or expand the chain of characters for the sequence*/
                tS[nbSeq-1].seq=(char*)realloc(tS[nbSeq-1].seq, sizeof(char)*(tS[nbSeq-1].lg+lgBuff+1));
                if(tS[nbSeq-1].seq==NULL){
                    fprintf(stderr, "readFasta:: Pb Allocation seq\n");
                    exit(1);
                }
                tS[nbSeq-1].seq[tS[nbSeq-1].lg]='\0';
                strcat(tS[nbSeq-1].seq, ligne);
                tS[nbSeq-1].lg+=lgBuff;
            }
        }
    }
    fclose(pFi);
    *pNbSeq=nbSeq;
    return tS;
}

/*-----------------------------------------------------------------------------------------------*/
/*											C(i,j)												 */
/*-----------------------------------------------------------------------------------------------*/

int** remplir_C(tgMotif unMotif, int seq_a_exclure, tgSeq *lesSeq, int numSeq) 
{
//This function constructs the C(i,j) table, which is the number of residue j at position
	//i in the motif. In the case of our example iteration, C(i,j) will look like this:

	/*

		C(i,j)				
	A	0.0	0.0	1.0	0.0
	B	1.0	1.0	0.0	1.0
	C	1.0	1.0	1.0	0.0
	D	0.0	0.0	0.0	1.0

	*/
	
	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. int seq_a_exclure, the number of a sequence that will be excluded from inclusion in the table. 
		//3. tgSeq *lesSeq, a table of tgSeq types that is our table of sequences
		//4. int numSeq, the number of sequences that we're dealing with 

	//Outputs:
		//1. table C(i,j)

	//créer le C(i,j) tableau. Ces dimensions sont [nombre de residues][longeur du motif cible].
	//Create the C(i,j) table. Its dimensions are [number of residues][length of the target motif]
	int **C = malloc(sizeof(int*) * NUM_RESIDUES);

	for(int i = 0; i < NUM_RESIDUES; i++) {
		C[i] = malloc(sizeof(int) * unMotif.lgMotif);
	}

	//remplir C(i,j) avec 0's. 
	for(int i=0; i < NUM_RESIDUES; i++) {
		for(int j=0; j < unMotif.lgMotif; j++) {
			C[i][j] = 0;
		}
	}

	//Ici on a une double boucle qui réitère au-dessus de chaque séquence (sauf l'un qui n'est pas inclus)
		//et puis elle réitère au-dessus de chaque résidue du motif, en réitérant les valeurs du tableau 
		//C(i,j).
	//Here we have a double loop which iterates over each sequence (except for the excluded one) and 
		//then it iterats over each residue of the motif, establishing the values of table C(i,j)
	for(int i=0; i < numSeq; i++) {
		for(int j=0; j < unMotif.lgMotif; j++) {
			//Cet "if-statement" exclut quelque séquence qui est choisit par seq_a_exclure
			//This if-statement excludes a sequence which is chosen by seq_a_exclure
			if(i != seq_a_exclure - 1) {
				//res égal à un des characters dans les motifs actuels. 
				//res is one of the characters in the current motifs
				char res = lesSeq[i].seq[unMotif.tDebut[i]+j];
				//calculer le "score" de la résidue, entre 1 et 26.
				//calculate the "score" of the residue, between 1 and 26 (regardless of AA or DNA)
				int indeceLettre = res - 'A';
				//itérer le valuer de C au rang qui correspond au résidue et au position actuel. 
				//Iterate the value in the row which corresponds to the residue and in the current position.
				C[indeceLettre][j]++;

			}
		}
	}

	return C;
}

/*-----------------------------------------------------------------------------------------------*/
/*											b(i)												 */
/*-----------------------------------------------------------------------------------------------*/

int* remplir_b(tgMotif unMotif, int seq_a_exclure, tgSeq *lesSeq, int numSeq) 
{		
//table b(i) is the number of each residue i which is observed outside of the motifs in each sequence. 
	//In the case of our example iteration, b(i) will look like this:

	/*

	   b(i)
	A  2.0
	B  3.0
	C  1.0
	D  3.0

	*/

	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. int seq_a_exclure, the number of a sequence that will be excluded from inclusion in the table. 
		//3. tgSeq *lesSeq, a table of tgSeq types that is our table of sequences
		//4. int numSeq, the number of sequences that we're dealing with 

	//Outputs:
		//1. table b, the counts of each residue outside of the motifs in each sequence. 


	//créer le b(i) tableau. La dimension est [nombre de residues].
	//Create the b(i) table. The length is [number of residues]
	int *b = malloc(sizeof(int*) * NUM_RESIDUES);

	//Puis, remplir b(i) avec 0's.
	//Fill b(i) with 0's.  
	for(int i=0; i < NUM_RESIDUES; i++) {
			b[i] = 0;
	}

	//Ici on a une double boucle qui réitère de chaque séquence (sauf l'un qui n'est pas inclus).
		//et puis elle réitère en dehors de chaque motif, en réitérant les valeurs du tableu b(i).
	//Here we want a double loop which iterates over each sequence (except for the sequence excluded),
		//and then it iterates over all positions outside of the motif in each sequence, counting up
		//the occurrences of each type of residue. 
	for(int i=0; i < numSeq; i++) {
		for(int j=0; j < lesSeq[i].lg; j++) {
			//Cet "if-statement" exclut quelque séquence qui est choisit par seq_a_exclure
			//This if-statement excludes the sequence chosen by seq_a_exclure.
			if(i != seq_a_exclure - 1) {
				//Puisqu'on cherche les résidues en dehors des motifs, il faut éviter les résidues
					//à l'intérieur des motifs
				//Because we're looking only at residues outside of the motifs, we have to avoid the residus
					//inside the motifs.
				if(unMotif.tDebut[i] > j || unMotif.tDebut[i]+unMotif.lgMotif <= j) {
					char res = lesSeq[i].seq[j];
					//calculer le "score" de la résidue, entre 1 et 26.
					//calculate the "score" of the residue (regardless of AA vs DNA), between 1 and 26. 
					int indeceLettre = res - 'A';
					//itérer le valuer de b au rang qui correspond au résidue. 
					//Iterate the value in table b at the row which corresponds to the residue observed. 
					b[indeceLettre]++;
				}
			}
		}
	}
	return b;
}

/*-----------------------------------------------------------------------------------------------*/
/*											rho(i)												 */
/*-----------------------------------------------------------------------------------------------*/

float* remplir_rho(tgMotif unMotif, tgSeq *lesSeq, int numSeq) 
{
//Table rho(i) is the ratio of each each residue i observed in all of the sequences (including the
	//excluded sequence and the motifs). In the case of our example iteration, rho(i) will look 
	//like this:

	/*

	   rho(i)
	A  0.2
	B  0.32
	C  0.2
	D  0.28

	*/

	// Inputs:
	// 		1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
	// 				along with other information.
	// 		2. tgSeq *lesSeq, a table of tgSeq types that is our table of sequences
	// 		3. int numSeq, the number of sequences that we're dealing with 

	// 	Outputs:
	// 		1. table rho(i), a table of the ratios of each residue i observed in all of the sequences

	//créer le rho(i) tableau. La dimension est [nombre de residues].
	//Create the rho(i) table. The dimension is [number of residues]
	float *rho = malloc(sizeof(float*) * NUM_RESIDUES);

	//Puis, remplir rho(i) avec 0's. 
	//Then, fill rho(i) with 0's.
	for(int i=0; i < NUM_RESIDUES; i++) {
			rho[i] = 0;
	}

	//Ici on a une double boucle qui réitère de chaque séquence (avec l'un qui n'est pas inclus).
		//et puis elle réitère sur chaque séquence, en réitérant les valeurs du tableu rho(i).
	//Here we have a double-loop which iterates over all sequences (including the excluded sequence)
		//then it iterates over each sequence and counts the occurrences of each residue. 
	for(int i=0; i < numSeq; i++) {
		for(int j=0; j < lesSeq[i].lg; j++) {
			//On utilise toutes les séquences données, donc il n'y a pas de if-statement ici.
			//We use all of the given sequences, so there's no if-statement here like the
				//previous functions. 
			char res = lesSeq[i].seq[j];
			//calculer le "score" de la résidue, entre 1 et 26.
			//calculate the "score" of the current residue, between 1 and 26 (regardless of AA or DNA)
			int indeceLettre = res - 'A';
			//itérer le valeur de rho au rang qui correspond au résidue. 
			//Iterate the value in rho(i) at the row which corresponds to the curren residue. 
			rho[indeceLettre]++;
		}
	}

	//Trouver la somme de la table rho(i):
	//Find the sum of the table rho(i)
	int somme = 0;
	for(int i=0; i < NUM_RESIDUES; i++) {
		somme += rho[i];
	}

	//Ici on a une boucle qui réitère au-dessus de la table rho(i) et calculer le ratio de chaque
		//résidue par raport à la somme. 
	//Here we have a loop which re-iterates over the table rho() to calcualte the ratio of each 
		//residue relative to the sum. 
	for(int i=0; i < NUM_RESIDUES; i++) {
		rho[i] = rho[i]/somme;
	}
	return rho;
}

/*-----------------------------------------------------------------------------------------------*/
/*											q(i,j)												 */
/*-----------------------------------------------------------------------------------------------*/

float** remplir_q(tgMotif unMotif, int **C, float *rho, int numSeq) 
{
//Table q(i,j) holds the frequency of residue j in position i of the motif, plus the "pseudocount".
	//See Lawrence et al. 1993 for further explanation of this table. In the case of our example
	//iteration, q(i,j) might look like this:

	// 		q(i,j)				
	// A	0.092820323	0.092820323	0.360769515	0.092820323
	// B	0.416461709	0.416461709	0.148512517	0.416461709
	// C	0.360769515	0.360769515	0.360769515	0.092820323
	// D	0.129948452	0.129948452	0.129948452	0.397897645

	//Note that the columns of q(i,j) sum to 1.

	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. int **C, created from the remplir_C() function above
		//3. int numSeq, the number of sequences that we're dealing with 

	//Outputs:
		//1. table q(i,j)
	
	//créer le q(i,j) tableau. Ces dimensions sont [nombre de residues][longeur du motif cible].
	//Create the q(i,j) table. The dimensions are [number of residues][length of the target motif]
	float **q = malloc(sizeof(float*) * NUM_RESIDUES);

	for(int i = 0; i < NUM_RESIDUES; i++) {
		q[i] = malloc(sizeof(float) * unMotif.lgMotif);
	}

	//Le valeur de beta est utilisé dans la formule si-dessous.
	//The value of beta is used in the formula below. 
	float beta = sqrt(numSeq);

	//utiliser la formule q(i,j) = C(i,j) + rho(i) * beta
		//						   --------------------- 
		//							    N - 1 + beta 
		//pour remplir q(i,j)
	//Use the above formula to fill q(i,j). 
	for(int i=0; i < NUM_RESIDUES; i++) {
		for(int j=0; j < unMotif.lgMotif; j++) {
			q[i][j] = (C[i][j] + 
					(rho[i] *
					beta)) /
					(numSeq -
					1 +
					beta);
		}
	}

 	return q;
}

/*-----------------------------------------------------------------------------------------------*/
/*											p(i)    											 */
/*-----------------------------------------------------------------------------------------------*/

float* remplir_p(tgMotif unMotif, tgSeq *lesSeq, int *b, float *rho, int numSeq, int seq_a_exclure) 
{
//p(i) is a one-dimensional analog of q(i,j). Given the use of b(i) in the formula
	//for each value of p(i), one can see that each value of p(i) is the general probability of 
	//a particular residue occurring at any position outside the motif. In the case of our example
	//iteration, p(i)looks like the following:

	/*

		p(i)
	A	0.218635767
	B	0.33118146
	C	0.125456931
	D	0.324725841

	*/

	//Note that p(i) sums to 1. 

	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. tgSeq *lesSeq, the table of sequences inputted
		//3. int *b, the table created by the function remplir_b()
		//4. float *rho, the table created by the function remplir_rho()
		//5. int numSeq, the number of sequences that we're dealing with 
		//6. int seq_a_exclure, the index in lesSeq of the sequence we want to exclude

	//Outputs:
		//1. table p(i)

	//créer le p(i) tableau. La dimension est [nombre de residues].
	//Create the p(i) table. The dimension is [number of residues].
	float *p = malloc(sizeof(float*) * NUM_RESIDUES);

	//utiliser la formule p(i) = b(i) + rho(i) * beta
		//						 --------------------- 
		//							 N - 1 + beta 
		//pour remplir p(i,j). C'est important de se souvenir, pourtant, que N soit différent pour 
		//p(i) que pour q(i,j). Ici, N est égal à la nombre de positions en dehors du motif dans les
		//séquences qu'on inclus. 
	//Use the above formula to calculate each value in p(i). It is important to remember, however,
		//that N is different for p(i) than it is foir q(i,j). Here, N is equal to the number of
		//position outside of the motifs in each of the sequences we've included. 

	//Définir N:
	//Define N:
	int N = 0;
	for(int i = 0; i < numSeq; i++) {
		if(i != seq_a_exclure - 1) {
			N += (lesSeq[i].lg - unMotif.lgMotif);
		}
	}

	float beta = sqrt(numSeq);

	for(int i=0; i < NUM_RESIDUES; i++) {
			p[i] = (b[i] + 
					(rho[i] *
					beta)) /
					(N +
					beta);
	}

	return p;
}

/*-----------------------------------------------------------------------------------------------*/
/*											F    												 */
/*-----------------------------------------------------------------------------------------------*/

float calculer_F(tgMotif unMotif, int **C, float **q, float *p) 
{
//The F score is the score of the alignment, calculated according to the formula outlined below,
	//from Lawrence et al. 1993. 

	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. int **C, the table created by the remplir_C() function
		//3. float **q, the table created by the remplir_q() function
		//4. float *p, the table created by the remplir_p() function. 

	//Outputs:
		//The score F of the alignment. 

	//Il faut calculer le "score" de l'alignement actuelle, par utiliser la formule suivante:
	//We have to calculate the score of the current alignment, by using the following formula. 

	/* 
			 w  NUM_RESIDUES
		    ___    ___
			\      \              q(i,j)
	F = 	/      /   C(i,j) log(----- )
			---    ---             p(j)
			i=1    j=1
	*/
	float F = 0;

	for(int i = 0; i < unMotif.lgMotif; i++) {
		for(int j = 0; j < NUM_RESIDUES; j++) {
			//beaucoup de valuers dans les tableaux q(i,j) et p(i) sont 0, donc il faut 
				//eviter log(0) et division par 0;
			//Many values in the q(i,j) and p(i) tables are 0, so we have to avoid 
				//log(0) and division by 0:
			if(q[j][i] >= .00001 && p[j] >= .00001) {
				F += (float)C[j][i] * log(q[j][i]/p[j]);
			}
		}
	}

	return F;
}

/*-----------------------------------------------------------------------------------------------*/
/*											Q(x)    											 */
/*-----------------------------------------------------------------------------------------------*/

float* remplir_Q(int nMotif, tgMotif unMotif, tgSeq *lesSeq, float **q, int seq_a_exclure)
{
//Q(x) holds the probabilities of generating each motif x according to the current pattern 
	//probabilities q(i,j). In the case of our example iteration, Q(x) looks like this:

	/*

		A	    D	    B	    C	    B	    D	 D	  A
	Qx	1.7E-04	8.1E-03	8.9E-03	7.8E-03	6.5E-04	0    0    0

	*/		 

	//Inputs:
		//1. nMotif, the number of possible positions in the target sequence 
		//2. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//3. tgSeq *lesSeq, the table of sequences inputted
		//4. float **q, the table created by the remplir_q() function
		//5. int seq_a_exclure, the index in lesSeq of the sequence we want to exclude

	//Outputs:
		//Table Q(x)

	//créer un tableau Q(x) pour porter la probabilité de chaque segment x selon les probabilités de
		//pattern dans q(i,j). La dimension va être Q[nombre de motifs possible dans la séquence cible]
	//Create a table Q(x) to hold the probability of each motif x accordning to the current pattern 
		//probabilities in q(i,j). The dimension is the number of motifs possible in the target sequence. 
	float *Q = malloc(sizeof(float*) * nMotif);

	for(int i = 0; i < nMotif; i++) {
		//remplir Q[x] avec le produit de toute les probabilités porté par q(i,j), par réitèrer
			//au-dessus de chaque résidue dans le motif actuel et référencer la résidue dans 
			//cette position dans q(i,j)
		//Fill Q(x) with the product of all of the probabilities held in q(i,j), by iterating
			//over each residue in the current motif and referencing the residue in that position
			//in q(i,j)
		Q[i] = 1;
		for(int j = 0; j < unMotif.lgMotif; j++) {
			char res = lesSeq[seq_a_exclure-1].seq[i+j];

			//calculer le "score" de la résidue, entre 1 et 26.
			//Calculate the "score" of the residue, between 1 and 26 (regardless of AA or DNA)
			int indeceLettre = res - 'A';

			//multiplié la valeur de Q[i] par la valeur correcte dans q[i,j]
			//Multiply the value of Q(i) by the correct value in q(i,j)
			Q[i] *= q[indeceLettre][j];
		}
	}

	return Q;
}

/*-----------------------------------------------------------------------------------------------*/
/*											P(x)    											 */
/*-----------------------------------------------------------------------------------------------*/

float* remplir_P(int nMotif, tgMotif unMotif, tgSeq *lesSeq, float *p, int seq_a_exclure)
{	
//P(x) holds the probabilities of generating each motif x according to the background probabilities
	//in p(i). In the case of our example iteration, P(x) looks like this:

	/*

		A	    D	    B	    C	    B	    D	 D	  A
	Px	2.9E-03	4.5E-03	4.5E-03	4.4E-03	7.6E-03	0    0    0

	*/

	//Inputs:
		//1. nMotif, the number of possible positions in the target sequence 
		//2. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//3. tgSeq *lesSeq, the table of sequences inputted
		//4. float *p, the table created by the remplir_p() function
		//5. int seq_a_exclure, the index in lesSeq of the sequence we want to exclude

	//Outputs:
		//Table P(x)

	//créer un tableau P[x] pour porter la probabilité de chaque segment x selon les probabilités de
		//pattern dans p[i]. La dimension va être P[nombre de motifs possible dans la séquence cible]
	//Create a table P(x) to store the probability of each motif x according to the pattern probabilities 
		//in p(i). The dimension is the number of motifs possible in the target sequence. 
	float *P = malloc(sizeof(float*) * nMotif);


	for(int i = 0; i < nMotif; i++) {
		//remplir P[x] avec le produit de toute les probabilités porté par p(i), par réitèrer
			//au-dessus de chaque résidue dans le motif actuel et référence la résidue dans p(i)
		//Fill P(x) with the product of all of the probabilities held in p(i), by iterating
			//over each residue in the current motif and referencing that residue in p(i)
		P[i] = 1;
		for(int j = 0; j < unMotif.lgMotif; j++) {
			char res = lesSeq[seq_a_exclure-1].seq[i+j];

			//calculer le "score" de la résidue, entre 1 et 26.
			//Calculate the "score" of the residue, between 1 and 26 (regardless of AA or DNA)
			int indeceLettre = res - 'A';

			//multiplié la valeur de P[i] par la valeur correcte dans p[i]
			//Multiply the value in P(i) by the correct value in p(i)
			P[i] *= p[indeceLettre];
		}
	}

	return P;
}

/*-----------------------------------------------------------------------------------------------*/
/*											A(x)    											 */
/*-----------------------------------------------------------------------------------------------*/

float* remplir_A(float *Q, float *P, int nMotif)
{	
//Table A(x) is simply the weight, or the probablity given to, each possible position for the 
	//motif in the target sequence. In the case of our example iteration, A(x) looks like this:

	/*

		A	    D	    B	    C	    B	    D	 D	  A
	Ax	0.056	1.820	1.987	1.773	0.085   0    0    0

	*/

	//Inputs:
		//1. float *Q, the table created in remplir_Q()
		//2. float *P, the table created in remplir_P()
		//3. int nMotif, the number of possible positions for the motif in the target sequence. 

	//Outputs:
		//1. Table A(x)

	float *A = malloc(sizeof(float*) * nMotif);

	//On veux aussi normaliser A(x), donc on calcule la somme aussi:
	//We also want to normalise A(x), so we also calculate the sum:
	float sommeA = 0;

	for(int i = 0; i < nMotif; i++) {
		A[i] = Q[i]/P[i];
		sommeA += A[i];
	}

	//Maintenant on normalise:
	//Now we normalise:
	for(int i = 0; i < nMotif; i++) {
		A[i] /= sommeA;
	}

	//Maintenant on change les valeurs de A(x) pour qu'on puisse choisir les indices
		//aléatoirement
	//Now we change the values of A(x) so that we can chose an index randomly 
		//using a random number, according to the probability distribution of A(x)
	for(int i = 1; i < nMotif; i++) {
		A[i] += A[i-1];
	}

	return A;
}

/*-----------------------------------------------------------------------------------------------*/
/*				calculer_F_pour_phaseShift (calculate F for the phase shift) 					 */
/*-----------------------------------------------------------------------------------------------*/

float calculer_F_pour_phaseShift(tgMotif unMotif, tgSeq *lesSeq, int numSeq, int shiftFactor) 
{
//This function is a helper function for the phaseShift() function below. It calculates the F score
	//of any alignment shifted from the current alignment by a value of shiftFactor. 

	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. tgSeq *lesSeq, the table of sequences inputted
		//3. int numSeq, the sumber of sequences inputted
		//4. int shiftFactor, the number of positions by which the current alignment should be shifted
				//in order to calculate the F score. 

	//Outputs:
		//1. F, the score of shifted alignment. 
		
	//itérer au-dessus de chaque tDebut dans unMotif et change la valeur par le shiftFactor
	//Iterate over each tDebut in unMotif and change its value by the value of shiftFactor

	for(int seq = 0; seq < numSeq; seq++) {

		unMotif.tDebut[seq] += shiftFactor;

		//si le tDebut est à une position qui n'est pas possible, la probabilité de ce shift 
			//est automatiquement 0, donc on retourne 0. 
		//If the tDebut is a position that is not possible, the probability of this shift is 
			//automatically 0, so we simply return 0. 
		if(unMotif.tDebut[seq] < 0 || unMotif.tDebut[seq] > lesSeq[seq].lg - unMotif.lgMotif) {
			return 0;
		}
	}

	/*REMPLIR C*/
	//on ne veut pas exclure une séquence, donc le "seq_a_exclure" est un seq qui n'existe pas
	//We don't want to exclude a sequence, so "seq_a_exclure" is a sequence that doesn't exist.
	int** C = remplir_C(unMotif, numSeq+1, lesSeq, numSeq);

	/*REMPLIR b*/
	int* b = remplir_b(unMotif, numSeq+1, lesSeq, numSeq);

	/*REMPLIR rho*/
	float* rho = remplir_rho(unMotif, lesSeq, numSeq);

	/*REMPLIR q*/
	float** q = remplir_q(unMotif, C, rho, numSeq);

	/*REMPLIRE p*/
	float* p = remplir_p(unMotif, lesSeq, b, rho, numSeq, numSeq+1);

	/*calculer F*/
	float F = calculer_F(unMotif, C, q, p);

	free(C);
	free(b);
	free(rho);
	free(q);
	free(p);

	return F;
}

/*-----------------------------------------------------------------------------------------------*/
/*											phaseShift    										 */
/*-----------------------------------------------------------------------------------------------*/

float* phaseShift(tgMotif unMotif, tgSeq *lesSeq, int maxShift, int numSeq)
{
//The purpose of phase shifting is described in Lawrence et al. 1993 as a method of escaping local
	//maxima while maintaining the general shape of the alignment. The idea is that the best
	//alignment could be the current alignment shifted right or left by one or several positions, 
	//and this is a great way to find it. 

	//Inputs:
		//1. tgMotif unMotif, a type that holds the starting positions of the motif in each sequence,
				//along with other information.
		//2. tgSeq *lesSeq, the table of sequences inputted
		//3. int maxShift, the maximum distance left or right that the alignment can be shifted.
		//4. int numSeq, the number of sequences inputted 

	//Outputs:
		//1. A table of two values:
			//1. the value of the phase shift chosen between -maxShift and maxShift. 
			//2. the F-score of the new alignment after the phase shift

	//Le tableau adjacent_Fs a un longueur de 2 fois maxShift, et il va guarder les scores 
		//F de chaque alignement adjacent à l'alignement actuel si on décale tous les débuts
		//par le même nombre de positions et dans la même direction. 
	//The table adjacent_Fs has a length of 2*maxShift, and it will store the F-score of each 
		//alignment adjacent to the current alignment if we shift all of the start positions
		//by the same number of position and in the same direction. 
	float *adjacent_Fs = malloc(sizeof(float*) * (2*maxShift));

	//guarder les tDebuts actuels dans un tableau pourqu'on puisse continuer à les utiliser
		//comme référence pour le shift.
	//Keep the current tDebuts in the table so that we can continue to use them as a reference
		//for the shift. 
	int *tDebuts_actuels = malloc(sizeof(int*) * numSeq);
	for(int i = 0; i < numSeq; i++) {
		tDebuts_actuels[i] = unMotif.tDebut[i];
	}

	//ce float est pour calculer la somme du tableau adjacent_Fs pour qu'on puisse le normaliser plus tard.
	//This float is to calculate the sum of the adjacent_Fs table so we can normalise later.  
	float somme_adjacent_Fs = 0;

	//remplir le tableau adjacent_Fs avec le F score de tous les alignements adjacents de 
		//l'alignement actuel. Après, on va normaliser la tableau pour choisir un phase shift 
		//selon la distribution des F scores.
	//Fill the table adjacent_Fs with the the F score of all of the alignments adjacent to the 
		//current alignment. Afterwards, we'll normalise the table in order to chose a 
		//phase shift according to the distribution of adjacent_Fs.  
	for(int i = -(maxShift); i < maxShift; i++) {

		if(i==0) {
			adjacent_Fs[i] = 0;
		} else {
			adjacent_Fs[i+maxShift] = calculer_F_pour_phaseShift(unMotif, lesSeq, numSeq, i);
			somme_adjacent_Fs += adjacent_Fs[i+maxShift];

			//Étant donné que calculer_F_pour_phaseShift change les tDebuts dans le struct unMotif, 
				//il faut les restorer par utiliser les valeurs guardés dans le tableau tDebuts_actuels. 
			//Given that calculer_F_pour_phaseShift changes the tDebuts in the struct unMotif, 
				//we need to restore them by using the values stored in the table tDebuts_actuels. 
			for(int j = 0; j < numSeq; j++) {
				unMotif.tDebut[j] = tDebuts_actuels[j];
			}
		}
	}

	//créer un autre tableau pour normaliser adjacent_Fs parce qu'on veut guarder les valeurs 
		//de F pour retourner de la fonctionne. 
	//Create another table to normalise adjacent_Fs because we want to keep the values of F
		//to return from this function. 
	float *adjacent_Fs_normalise = malloc(sizeof(float*) * (2*maxShift));

	//normaliser le tableau adjacent_Fs
	//normalise the table adjacent_Fs
	for(int i = 0; i < (2*maxShift); i++) {
		adjacent_Fs_normalise[i] = adjacent_Fs[i]/somme_adjacent_Fs;
	}
	for(int i = 1; i < (2*maxShift); i++) {
		adjacent_Fs_normalise[i] += adjacent_Fs_normalise[i-1];
	}

	//générer une numéro aléatoire pour choisir une des indices de adjacent_Fs_normalise.
	//Generate a random number to chose one of the indices in adjacent_Fs_normalise.
	double aleatoire_pour_F = rand()/(double)(RAND_MAX);

	float shiftChoisi = 0;
	float fChoisi = 0;
	//choisir une des indices de adjacent_Fs_normalise basée sur la numéro aléatoire.
	//Chose one of the indices of adjacent_Fs_normalise based on the probability 
		//distribution of adjacent_Fs_normilise, basedo on the random number. 
	for(int i = -(maxShift); i < maxShift; i++) {
		if(adjacent_Fs_normalise[i+maxShift] > aleatoire_pour_F) {
			shiftChoisi = i;
			fChoisi = adjacent_Fs[i+maxShift];
			break;
		}
	}

	float* returnTableau = malloc(sizeof(float*) * (2));
	returnTableau[0] = shiftChoisi;
	returnTableau[1] = fChoisi;

	free(adjacent_Fs);
	free(tDebuts_actuels);
	free(adjacent_Fs_normalise);

	return returnTableau;
}




