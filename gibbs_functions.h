#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>


//D'abord, on crée des nouveaux types, qui s'appellent tgSeq et tgMotif.

//tgSeq est une collection d'information sur une seule séquence, qu'on peut
	//utiliser pour créer un tableau de séquences
typedef struct _tgSeq { char *titre ; char *seq ; int lg ; } tgSeq ;

//tgMotif est une collection d'information sur une seule motif, qu'on peut manipuler 
	//pour éventuellement trouver le cadre séquence. 
typedef struct _tgMotif { int lgMotif ; int nbSeq ; tgSeq *tSeq ; int *tDebut ; } tgMotif;

//Ici, on défine le nombre de résidues. C'est important car on peut utiliser cet algorithme 
	//pour les séquences de protéines et aussi de ADN.
#define NUM_RESIDUES 26
#define NUM_ITERATIONS 1000
#define NUM_COMMANDS 2
#define MAX_SEQ_LN 1000
#define LG_MOTIF 6
#define PHASE_SHIFT_FREQUENCY 1001
#define MAX_SHIFT 5


/*-----------------------------------------------------------------------------------------------*/
/*											readFasta											 */
/*-----------------------------------------------------------------------------------------------*/

tgSeq* readFasta(char *nomFi, int *pNbSeq)
{

    tgSeq *tS=NULL;

    char ligne[MAX_SEQ_LN];
    int lgBuff;
    int nbSeq=0;
    FILE *pFi;

    pFi=fopen(nomFi, "rt"); /*Ouvrir le fichier */
    if(pFi==NULL){
        fprintf(stderr, "readFasta:: fichier introuvable: %s\n", nomFi);
        return NULL;
    }

    /*Lire les  lignes*/
    while(fgets(ligne, MAX_SEQ_LN, pFi)!=NULL) {
    	/*Calculer la longueur de la ligne*/
        lgBuff=strlen(ligne); 
        /*Enlever le retour à la ligne s'il y en a un*/
        if(ligne[lgBuff-1]=='\n'){ 
            ligne[lgBuff-1]='\0';
            lgBuff--;
        }
        /*Si la ligne n'est pas vide*/
        if (lgBuff>0) { 
            /*Si c'est une ligne de titre*/
            if(ligne[0]=='>') { 
                nbSeq++;               
                /*Ajouter une case au tableau des séquences*/
                tS=(tgSeq*)realloc(tS, sizeof(tgSeq)*(nbSeq));               
                if(tS==NULL){
                    fprintf(stderr, "readFasta:: Erreur Allocation tS\n");
                    exit(1);
                }
                /*créer les membres de la nouvelle séquence*/
                tS[nbSeq-1].lg=0;
                tS[nbSeq-1].seq=NULL;
                /*Allouer une chaine de caracteres pour le titre */
                tS[nbSeq-1].titre=(char*)malloc(sizeof(char)*(lgBuff+1)); 
                if(tS[nbSeq-1].titre==NULL){
                    fprintf(stderr, "readFasta:: Pb Allocation titre\n");
                    exit(1);
                }
                strcpy(tS[nbSeq-1].titre, ligne); /*Recopier la ligne de titre*/
            /*C'est une ligne de séquence*/
            } else {
                if(nbSeq==0){
                    fprintf(stderr, "readFasta:: format incorrect, exiting\n");
                    exit(1);
                }
                /*Allouer ou agrandir la chaine de caracteres pour la sequence */
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

	//créer le C(i,j) tableau. Ces dimensions sont [nombre de residues][longeur du motif cible].
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
	for(int i=0; i < numSeq; i++) {
		for(int j=0; j < unMotif.lgMotif; j++) {
			//Cet "if-statement" exclut quelque séquence qui est choisit par SEQ_A_EXCLURE
			if(i != seq_a_exclure - 1) {
				//res égal à un des characters dans les motifs actuels. 
				char res = lesSeq[i].seq[unMotif.tDebut[i]+j];
				//calculer le "score" de la résidue, entre 1 et 26.
				int indeceLettre = res - 'A';
				//itérer le valuer de C au rang qui correspond au résidue et au position actuel. 
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
	
	//créer le b(i) tableau. La dimension est [nombre de residues].
	int *b = malloc(sizeof(int*) * NUM_RESIDUES);

	//Puis, remplir b(i) avec 0's. 
	for(int i=0; i < NUM_RESIDUES; i++) {
			b[i] = 0;
	}

	//Ici on a une double boucle qui réitère de chaque séquence (sauf l'un qui n'est pas inclus).
		//et puis elle réitère en dehors de chaque motif, en réitérant les valeurs du tableu C(i).
	for(int i=0; i < numSeq; i++) {
		for(int j=0; j < lesSeq[i].lg; j++) {
			//Cet "if-statement" exclut quelque séquence qui est choisit par SEQ_A_EXCLURE
			if(i != seq_a_exclure - 1) {
				//Puisqu'on cherche les résidues en dehors des motifs, il faut éviter les résidues
					//à l'intérieur des motifs
				if(unMotif.tDebut[i] > j || unMotif.tDebut[i]+unMotif.lgMotif <= j) {
					char res = lesSeq[i].seq[j];

					//calculer le "score" de la résidue, entre 1 et 26.
					int indeceLettre = res - 'A';
					//itérer le valuer de b au rang qui correspond au résidue. 
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

float* remplir_rho(tgMotif unMotif, int seq_a_exclure, tgSeq *lesSeq, int numSeq) 
{

	//créer le rho(i) tableau. La dimension est [nombre de residues].
	float *rho = malloc(sizeof(float*) * NUM_RESIDUES);

	//Puis, remplir rho(i) avec 0's. 
	for(int i=0; i < NUM_RESIDUES; i++) {
			rho[i] = 0;
	}

	//Ici on a une double boucle qui réitère de chaque séquence (sauf l'un qui n'est pas inclus).
		//et puis elle réitère en dehors de chaque motif, en réitérant les valeurs du tableu rho(i).
	for(int i=0; i < numSeq; i++) {
		for(int j=0; j < lesSeq[i].lg; j++) {
			//On utilise toutes les séquences données, donc il n'y a pas de if-statement ici. 
			char res = lesSeq[i].seq[j];

			//calculer le "score" de la résidue, entre 1 et 26.
			int indeceLettre = res - 'A';
			//itérer le valuer de rho au rang qui correspond au résidue. 
			rho[indeceLettre]++;
		}
	}

	//Trouver la somme de la table rho(i):
	int somme = 0;
	for(int i=0; i < NUM_RESIDUES; i++) {
		somme += rho[i];
	}

	//Ici on a une boucle qui réitère au-dessus de la table rho(i) et calculer le ratio de chaque
		//rédidue par raport à la somme. 
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
	
	//créer le q(i,j) tableau. Ces dimensions sont [nombre de residues][longeur du motif cible].
	float **q = malloc(sizeof(float*) * NUM_RESIDUES);

	for(int i = 0; i < NUM_RESIDUES; i++) {
		q[i] = malloc(sizeof(float) * unMotif.lgMotif);
	}

	float beta = sqrt(numSeq);

	//utiliser la formule q(i,j) = C(i,j) + rho(i) * beta
		//						   --------------------- 
		//							    N - 1 + beta 
		//pour remplir q(i,j)
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

float* remplir_p(tgMotif unMotif, tgSeq *lesSeq, int *b, float *rho, int numSeq, int SEQ_A_EXCLURE) 
{
	//créer le p(i) tableau. La dimension est [nombre de residues].
	float *p = malloc(sizeof(float*) * NUM_RESIDUES);

	//utiliser la formule p(i) = b(i) + rho(i) * beta
		//						 --------------------- 
		//							 N - 1 + beta 
		//pour remplir p(i,j). C'est important de se souvenir, pourtant, que N soit différent pour 
		//p(i) que pour q(i,j). Ici, N est égal à la nombre de positions en dehors du motif dans les
		//séquences qu'on inclus. 
		//TODO: Créer une fonction externelle pour cette boucle. 

	//Définir N:
	int N = 0;
	for(int i = 0; i < numSeq; i++) {
		if(i != SEQ_A_EXCLURE - 1) {
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
	//ensuite, il faut calculer le "score" d'alignement actuelle, par utiliser la formule suivante:

	//float *F = malloc(sizeof(double));

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

float* remplir_Q(int nMotif, tgMotif unMotif, tgSeq *lesSeq, float **q, int SEQ_A_EXCLURE)
{
	//créer un tableau Q[x] pour porter la probabilité de chaque segment x selon les probabilités de
		//pattern dans q(i,j). La dimension va être Q[nombre de motifs possible dans la séquence cible]
	float *Q = malloc(sizeof(float*) * nMotif);

	for(int i = 0; i < nMotif; i++) {
		//remplir Q[x] avec le produit de toute les probabilités porté par q(i,j), par réitèrer
			//au-dessus de chaque résidue dans le motif actuel et référencer la résidue dans 
			//cette position dans q(i,j)
		Q[i] = 1;
		for(int j = 0; j < unMotif.lgMotif; j++) {
			char res = lesSeq[SEQ_A_EXCLURE-1].seq[i+j];

			//calculer le "score" de la résidue, entre 1 et 26.
			int indeceLettre = res - 'A';

			//multiplié la valeur de Q[i] par la valeur correcte dans q[i,j]
			Q[i] *= q[indeceLettre][j];
		}
	}

	return Q;
}

/*-----------------------------------------------------------------------------------------------*/
/*											P(x)    											 */
/*-----------------------------------------------------------------------------------------------*/

float* remplir_P(int nMotif, tgMotif unMotif, tgSeq *lesSeq, float *p, int SEQ_A_EXCLURE)
{		
	//créer un tableau P[x] pour porter la probabilité de chaque segment x selon les probabilités de
	//pattern dans p[i]. La dimension va être Q[nombre de motifs possible dans la séquence cible]
	float *P = malloc(sizeof(float*) * nMotif);


	for(int i = 0; i < nMotif; i++) {
		//remplir P[x] avec le produit de toute les probabilités porté par p(i), par réitèrer
			//au-dessus de chaque résidue dans le motif actuel et référence la résidue dans p(i)
		P[i] = 1;
		for(int j = 0; j < unMotif.lgMotif; j++) {
			char res = lesSeq[SEQ_A_EXCLURE-1].seq[i+j];

			//calculer le "score" de la résidue, entre 1 et 26.
			int indeceLettre = res - 'A';

			//multiplié la valeur de P[i] par la valeur correcte dans p[i]
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

	float *A = malloc(sizeof(float*) * nMotif);

	//On veux aussi normaliser A(x), donc on calcule la somme aussi:
	float sommeA = 0;

	for(int i = 0; i < nMotif; i++) {
		A[i] = Q[i]/P[i];
		sommeA += A[i];
	}

	//Maintenant on normalise:
	for(int i = 0; i < nMotif; i++) {
		A[i] /= sommeA;
	}

	//Maintenant on change les valeurs de A(x) pour qu'on puisse choisir les indices
		//aléatoirement

	for(int i = 1; i < nMotif; i++) {
		A[i] += A[i-1];
	}

	return A;

}

/*-----------------------------------------------------------------------------------------------*/
/*									calculer_F_pour_phaseShift    								 */
/*-----------------------------------------------------------------------------------------------*/

float calculer_F_pour_phaseShift(tgMotif unMotif, tgSeq *lesSeq, int numSeq, int shiftFactor) {

	//itérer au-dessus de chaque tDebut dans unMotif et change la valeur par le shiftFactor
	for(int seq = 0; seq < numSeq; seq++) {

		unMotif.tDebut[seq] += shiftFactor;

		//si le tDebut est à une position qui n'est pas possible, la probabilité de ce shift 
			//est automatiquement 0, donc on retourne 0. 
		if(unMotif.tDebut[seq] < 0 || unMotif.tDebut[seq] > lesSeq[seq].lg - unMotif.lgMotif) {
			return 0;
		}
	}

	/*REMPLIR C*/
	//on ne veut pas exclure une séquence, donc le "seq_a_exclure" est un seq qui n'existe pas
	int** C = remplir_C(unMotif, numSeq+1, lesSeq, numSeq);

	/*REMPLIR b*/
	int* b = remplir_b(unMotif, numSeq+1, lesSeq, numSeq);

	/*REMPLIR rho*/
	float* rho = remplir_rho(unMotif, numSeq+1, lesSeq, numSeq);

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
	
	//Le tableau adjacent_Fs a un longueur de 2 fois maxShift, et il va guarder les scores 
		//F de chaque alignement adjacent à l'alignement actuel si on décale tous les débuts
		//par le même nombre de positions et dans la même direction. 
	float *adjacent_Fs = malloc(sizeof(float*) * (2*maxShift));

	//guarder les tDebuts actuels dans un tableau pourqu'on puisse continuer à les utiliser
		//comme référence pour le shift.
	int *tDebuts_actuels = malloc(sizeof(int*) * numSeq);
	for(int i = 0; i < numSeq; i++) {
		tDebuts_actuels[i] = unMotif.tDebut[i];
	}

	//ce float est pour calculer la somme du tableau adjacent_Fs pour qu'on puisse le normaliser plus tard. 
	float somme_adjacent_Fs = 0;

	//remplir le tableau adjacent_Fs avec le F score de tous les alignements adjacents de 
		//l'alignement actuel. Après, on va normaliser la tableau pour choisir un phase shift 
		//selon la distribution des F scores. 
	for(int i = -(maxShift); i < maxShift; i++) {

		adjacent_Fs[i+maxShift] = calculer_F_pour_phaseShift(unMotif, lesSeq, numSeq, i);
		somme_adjacent_Fs += adjacent_Fs[i+maxShift];

		//Étant donné que calculer_F_pour_phaseShift change les tDebuts dans le struct unMotif, 
			//il faut les restorer par utiliser les valeurs guardés dans le tableau tDebuts_actuels. 
		for(int j = 0; j < numSeq; j++) {
			unMotif.tDebut[j] = tDebuts_actuels[j];
		}
	}

	//créer un autre tableau pour normaliser adjacent_Fs parce qu'on veut guarder les valeurs 
		//de F pour retourner de la fonctionne. 
	float *adjacent_Fs_normalise = malloc(sizeof(float*) * (2*maxShift));

	//normaliser le tableau adjacent_Fs
	for(int i = 0; i < (2*maxShift); i++) {
		adjacent_Fs_normalise[i] = adjacent_Fs[i]/somme_adjacent_Fs;
	}
	for(int i = 1; i < (2*maxShift); i++) {
		adjacent_Fs_normalise[i] += adjacent_Fs_normalise[i-1];
	}

	//générer une numéro aléatoire pour choisir une des indices de A(x).
	double aleatoire_pour_F = rand()/(double)(RAND_MAX);

	float shiftChoisi = 0;
	float fChoisi = 0;
	//choisir une indices de A(x) basée sur la numéro aléatoire.
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

	return returnTableau;

}




