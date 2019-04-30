#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "gibbs_functions.h"

/*-----------------------------------------------------------------------------------------------*/
/*											   											     	 */
/*											  MAIN  											 */
/*											    												 */
/*-----------------------------------------------------------------------------------------------*/

int main() 
{

	//Initialiser un clock pour calculer le runtime.
	clock_t start, end;
	double cpu_time_used;
    start = clock();

    //les lignes suivantes lisent la fichier fasta et créent un tableau de tgSeq avec des séquences
    int nbSeq = 0;
    tgSeq* lesSeq = readFasta("lipocalin.fst", &nbSeq);

	//initialiser un tgMotif qui s'appelle "unMotif", et le remplir avec les indices de debut qui
		//sont choisis "aléatoirement" (en ce cas, pas exactement aléatoirement)
	tgMotif unMotif;
	unMotif.lgMotif = LG_MOTIF;
	unMotif.nbSeq = nbSeq;
	unMotif.tSeq = lesSeq;
	unMotif.tDebut = malloc(sizeof(int)*unMotif.nbSeq);

	//important pour choisir les nombres aléatoires aléatoirement! 
	srand ( time(NULL) );

	//créer nombres aléatoires entre 0 et la longeur de chaque séquence, pour être les prémiers
		//débuts du motif dans chaque séquence. 
	for(int i = 0; i < nbSeq; i++) {
		int tDebut = rand() % (lesSeq[i].lg - unMotif.lgMotif + 1);
		unMotif.tDebut[i] = tDebut;
	}

	//afficher les débuts originaux
	printf("Originals: \n");
	printf("Seq1: %d, Seq2: %d, Seq3: %d \n", unMotif.tDebut[0], unMotif.tDebut[1], unMotif.tDebut[2]);

	printf("-------------------\n");

	//créer le variable de score, F_global, et le int pour le meilleur alignement:
	//F_global et meilleur_index suivent le meilleur alignement actuel, qui n'est pas nécessairement
		//le meilleur alignement trouvé. 
	float F_global = .0000001; 
	int meilleur_index = 0;

	//le tableau suivre_F guarde tous les valeurs de F_global et F pour qu'on puisse grapher les itérations. 
	double suivre_F[NUM_ITERATIONS][3];

	//meilleur_F et index_de_meilleur_F guardent le meilleur valeur de F jamais trouvé. 
	float meilleur_F = 0;
	int index_de_meilleur_F = meilleur_index;

	//le tableau "resultats" va porter les tDebut après chaque réitération
	int resultats[NUM_ITERATIONS][nbSeq+1];

	//remplir "resultats" avec les tDebuts actuels
	resultats[0][0] = 0;
	for(int j = 0; j < nbSeq; j++) { 
		resultats[0][j+1] = unMotif.tDebut[j];
	}

	/*maintenant on fait des itérations*/
	for(int iteration = 0; iteration < NUM_ITERATIONS; iteration++) 
	{
		
		if(iteration % PHASE_SHIFT_FREQUENCY == 0 && iteration > 0) {

			/*---------------------------------------------------------------------*/
			/*						OPTION 1: PHASE SHIFT 						   */
			/*---------------------------------------------------------------------*/
			/*Parfois, l'alignement devient coincé dans un maximum locale. Selon 
			  Lawrence et al. 1993, une solution possible pour ce problème est de 
			  faire un "phase shift", c'est-à-dire, décaler tous les débuts par le 
			  même nombre de positions. Le nombre de décalage est décidé par le
			  phaseShift fonctionne dans gibbs_functions.h, par calculer le score F
			  de chaque set de Ax adjacent de l'Ax actuel, au moins de MAX_SHIFT de 
			  l'Ax actuel. On fait un phaseShift chaque PHASE_SHIFT_FREQUENCY itérations*/

			float* shiftNombres = phaseShift(unMotif, lesSeq, MAX_SHIFT, nbSeq);

			//shiftBy est le nombre de décalage, et nouveauF est le score F du set de Ax 
				//choisi par la fonctionne phaseShift
			int shiftBy = (float)shiftNombres[0];
			float nouveauF = shiftNombres[1];

			//changer les debuts par le nombre shiftBy. 
			for(int i = 0; i < nbSeq; i++) {
				unMotif.tDebut[i] += shiftBy;
			}

			//mettre à jour le tableau resultats avec les tDebuts changés. 
			resultats[iteration][0] = iteration;
			for(int j = 0; j < nbSeq; j++) { 
				resultats[iteration][j+1] = unMotif.tDebut[j];
			}

			//guarder le vrai meilleur F. 
			if(nouveauF > meilleur_F) {
				meilleur_F = nouveauF;
				index_de_meilleur_F = iteration;
			}
			
			//Si on ne change pas meilleur_index ni F_global, l'algorithme reviendra au même maximum locale
				//qu'avant, donc il faux assurer que le nouveau F_global soit le nouveauF trouvé par phaseShift. 
			meilleur_index = iteration;
			F_global = nouveauF;

			//on veut aussi plotter F et F_global pour voir le performance de l'algorithme. On guarde les valeurs
			//de F et F_global dans le table suivant:
			suivre_F[iteration][0] = iteration;
			suivre_F[iteration][1] = nouveauF;
			suivre_F[iteration][2] = F_global;

		} else {

			/*---------------------------------------------------------------------*/
			/*						OPTION 2: GIBBS SAMPLING					   */
			/*---------------------------------------------------------------------*/
			/*Dans cette option, on emploi le processus de Gibbs Sampling d'après Lawrence et al. 1993.*/ 
			
			//définir un séquence aléatoirement à exlure. 
			int SEQ_A_EXCLURE = (rand() % nbSeq) + 1;

			/*-----------------------------------------------------------------------------------------------*/
			/*											C(i,j)												 */
			/*-----------------------------------------------------------------------------------------------*/

			int** C = remplir_C(unMotif, SEQ_A_EXCLURE, lesSeq, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											b(i)												 */
			/*-----------------------------------------------------------------------------------------------*/

			int* b = remplir_b(unMotif, SEQ_A_EXCLURE, lesSeq, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											rho(i)												 */
			/*-----------------------------------------------------------------------------------------------*/

			float* rho = remplir_rho(unMotif, SEQ_A_EXCLURE, lesSeq, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											q(i,j)												 */
			/*-----------------------------------------------------------------------------------------------*/

			float** q = remplir_q(unMotif, C, rho, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											p(i)    											 */
			/*-----------------------------------------------------------------------------------------------*/

			float* p = remplir_p(unMotif, lesSeq, b, rho, nbSeq, SEQ_A_EXCLURE);

			/*-----------------------------------------------------------------------------------------------*/
			/*											F    												 */
			/*-----------------------------------------------------------------------------------------------*/

			//calculer F
			float F = calculer_F(unMotif, C, q, p);

			//guarder le vrai meilleur F. 
			if(F > meilleur_F) {
				meilleur_F = F;
				index_de_meilleur_F = iteration;
			}
			
			/*---------------------*/
			/*CRITERE DE METROPOLIS*/
			/*---------------------*/
			//"probabilite" décroisse quand F << F_global.
			float probabilite = F/F_global;
			double metropolis_chiffre = rand()/(double)(RAND_MAX);
			
			resultats[iteration][0] = iteration;
			for(int j = 0; j < nbSeq; j++) { 
				resultats[iteration][j+1] = unMotif.tDebut[j];
			}

			//Avec la critère de Metropolis, le variable "meilleur_index" va rester la même sauf si le nouveau 
				//alignement a un mieux score que le meilleur alignement, ou si l'alignement actuel est choisi
				//aléatoirement, selon la probabilité calculé dans le variable "probabilite". 
			if(F > F_global || probabilite >= metropolis_chiffre) {
				F_global = F;
				meilleur_index = iteration;
			} 

			//Ici, on guarde l'alignement guardé dans le "meilleur_index" de resultats. Le meilleur_index peut 
				//être l'index de l'itération actuelle, ou l'index d'une itération passée qui a trouvé la 
				//meilleur alignement. 
			resultats[iteration][0] = iteration;
			for(int col = 0; col < nbSeq; col++) {
				unMotif.tDebut[col] = resultats[meilleur_index][col+1];
			}

			//remplir "resultats" avec les tDebuts actuels
			for(int j = 0; j < nbSeq; j++) { 
				resultats[iteration][j+1] = unMotif.tDebut[j];
			}

			//on veut aussi plotter F et F_global pour voir le performance de l'algorithme. On guarde les valeurs
			//de F et F_global dans le table suivant:
			suivre_F[iteration][0] = iteration;
			suivre_F[iteration][1] = F;
			suivre_F[iteration][2] = F_global;

			/*-----------------------------------------------------------------------------------------------*/
			/*											Q(x)    											 */
			/*-----------------------------------------------------------------------------------------------*/

			//nMotif réprésent le nombre de motifs possible de longeur unMotif.lgMotif qui sont dans
			//la séquence exclue. 
			int nMotif = lesSeq[SEQ_A_EXCLURE-1].lg - unMotif.lgMotif + 1;

			float* Q = remplir_Q(nMotif, unMotif, lesSeq, q, SEQ_A_EXCLURE);

			/*-----------------------------------------------------------------------------------------------*/
			/*											P(x)    											 */
			/*-----------------------------------------------------------------------------------------------*/

			float* P = remplir_P(nMotif, unMotif, lesSeq, p, SEQ_A_EXCLURE);

			//Maintenant on calcule les ratios entre Q(x) et P(x) pour créer A(x). 

			/*-----------------------------------------------------------------------------------------------*/
			/*											A(x)    											 */
			/*-----------------------------------------------------------------------------------------------*/
			
			float* A = remplir_A(Q, P, nMotif);

			
			//générer une numéro aléatoire pour choisir une des indices de A(x).
			double aleatoire = rand()/(double)(RAND_MAX);

			//printf("Aleatoire: %f \n", aleatoire);

			//choisir une indices de A(x) basée sur la numéro aléatoire.
			for(int i = 0; i < nMotif; i++) {
				if(A[i] > aleatoire) {
					unMotif.tDebut[SEQ_A_EXCLURE-1] = i;
					break;
				}
			}

			//vider la mémoire
			free(C);
			free(b);
			free(rho);
			free(q);
			free(p);
			free(Q);
			free(P);
			free(A);
		
		}
	}

	//on a guarder l'index de meilleur F, quelque part dans le tableau "resultats". Maintenant, les 
		//tDebuts de chaque séquence est finalement égal au tDebuts du meilleur alignement. 
	for(int col = 0; col < nbSeq; col++) {
		unMotif.tDebut[col] = resultats[index_de_meilleur_F][col+1];
	}
	
	/*-----------------------------------------------------------------------------------------------*/
	/*											AFFICHAGE											 */
	/*-----------------------------------------------------------------------------------------------*/

	//afficher toutes les informations importants. 
	printf("Meilleur_F: %f \n", meilleur_F);
	printf("F_global: %f \n", F_global);
	printf("index_de_meilleur_F: %d \n", index_de_meilleur_F);
	for(int i=0; i < nbSeq; i++) {
		int tDebut = resultats[index_de_meilleur_F][i+1];
		char *seq = lesSeq[i].seq;
		printf("Sequence %d motif result: %d \n", i, tDebut);
		printf("%.*s\n", unMotif.lgMotif, seq + tDebut);
	}

	/*-----------------------------------------------------------------------------------------------*/
	/*											PLOTTING											 */
	/*-----------------------------------------------------------------------------------------------*/
	
	char * commandsForGnuplot[] = {"set title \"Fglobal sur toutes les itérations\"", "plot 'data.temp' with lines"};
    FILE * temp = fopen("data.temp", "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    for(int i = 0; i < NUM_ITERATIONS; i++) {
		//fprintf(temp, "%-5f %-5f %-5f\n", suivre_F[i][0], suivre_F[i][1], suivre_F[i][2]);
		fprintf(temp, "%-5f %-5f\n", suivre_F[i][0], suivre_F[i][2]);
	}

    for (i=0; i < NUM_COMMANDS; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }

	//imprimer le runtime
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Runtime: %f seconds\n", cpu_time_used);

	return 1;

}