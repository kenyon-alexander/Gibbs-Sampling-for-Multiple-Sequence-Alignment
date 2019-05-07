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
	//Initialise the clock to calculate the runtime.
	clock_t start, end;
	double cpu_time_used;
    start = clock();

    //les lignes suivantes lisent la fichier fasta et créent un tableau de tgSeq avec des séquences
    //the following lines read the fasta file and create a table of tgSeq variables with the sequences
    	//loaded into them.
    int nbSeq = 0;
    tgSeq* lesSeq = readFasta("lipocalin.fst", &nbSeq);

	//initialiser un tgMotif qui s'appelle "unMotif", et le remplir avec les indices de debut qui
		//sont choisis "aléatoirement" (en ce cas, pas exactement aléatoirement)
    //Initialise a tgMotif variable called "unMotif", and fill it with the start indices chosen 
    	//randomly.
	tgMotif unMotif;
	unMotif.lgMotif = LG_MOTIF;
	unMotif.nbSeq = nbSeq;
	unMotif.tSeq = lesSeq;
	unMotif.tDebut = malloc(sizeof(int)*unMotif.nbSeq);

	//important pour choisir les nombres aléatoires aléatoirement!
	//We need to set the pseudo-random number generator's time to null so that the numbers can be chosen 
		//differently each time.
	srand ( time(NULL) );

	//créer nombres aléatoires entre 0 et la longeur de chaque séquence, pour être les prémiers
		//débuts du motif dans chaque séquence. 
	//Chose random numbers betwen 0 and the length of each sequence, to be the initial motif positions
		//in each sequence.
	for(int i = 0; i < nbSeq; i++) {
		int tDebut = rand() % (lesSeq[i].lg - unMotif.lgMotif + 1);
		unMotif.tDebut[i] = tDebut;
	}

	//afficher les débuts originaux
	//Print the original motif positions
	printf("Originals: \n");
	for(int i=0;i<nbSeq;i++) {
		printf("Seq%d: %d ", i+1, unMotif.tDebut[i]);
	}
	printf("\n");
	//printf("Seq1: %d, Seq2: %d, Seq3: %d \n", unMotif.tDebut[0], unMotif.tDebut[1], unMotif.tDebut[2]);

	printf("-------------------\n");

	//meilleur_F_actuel et meilleur_index suivent le meilleur alignement actuel, qui n'est pas nécessairement
		//le meilleur alignement trouvé. 
	//meilleur_F_actuel and meilleur_index track the current best alignment, which is not necessarily the best alignment 
		//found over all iterations. 
	float meilleur_F_actuel = .0000001; 
	int meilleur_index = 0;

	//meilleur_F et index_de_meilleur_F guardent le meilleur valeur de F jamais trouvé. 
	//meilleur_F and index_de_meilleur_F save the best score ever found for F. 
	float meilleur_F = 0;
	int index_de_meilleur_F = meilleur_index;

	//le tableau suivre_F guarde tous les valeurs de meilleur_F_actuel et F pour qu'on puisse grapher les itérations. 
	//the table suivre_F stores all of the values of meilleur_F_actuel and F so we can graph the progress of meilleur_F_actuel and F. 
	// static double suivre_F[NUM_ITERATIONS][3];
	double **suivre_F = malloc(sizeof(double*) * NUM_ITERATIONS);

	for(int i = 0; i < NUM_ITERATIONS; i++) {
		suivre_F[i] = malloc(sizeof(double) * 3);
	}

	//le tableau "resultats" va porter les tDebut après chaque réitération
	//The table "resultats" will hold the motif positions (tDebut) after each iteration. 
	// static int resultats[NUM_ITERATIONS][nbSeq+1];

	int **resultats = malloc(sizeof(int*) * NUM_ITERATIONS);

	for(int i = 0; i < NUM_ITERATIONS; i++) {
		resultats[i] = malloc(sizeof(int) * nbSeq+1);
	}

	//remplir "resultats" avec les tDebuts actuels
	//fill the first row of "resultats" with the current motif positions (tDebuts) 
	resultats[0][0] = 0;
	for(int j = 0; j < nbSeq; j++) { 
		resultats[0][j+1] = unMotif.tDebut[j];
	}

	/*maintenant on fait des itérations*/
	//Now we perform the iterations. 
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

			/*Sometimes, the alignment gets stuck in a local maximum. According to 
			  Lawrence et al. 1993, a possible solution for this problem is to perform
			  a "phase shift", in other words, to shift all of the motif positions (tDebuts)
			  by the same amount. The leght of the phase shift is decided by the
			  phaseShift() function in gibbs_functions.h, by calculating the F score 
			  of each set of tDebuts adjacent to the current set, up to the length of 
			  the value of the MAX_SHIFT constant. As it currently stands, the algorithm
			  performs a phase shift each PHASE_SHIFT_FREQUENCY interations*/ 

			
			float* shiftNombres = phaseShift(unMotif, lesSeq, MAX_SHIFT, nbSeq);

			//shiftBy est le nombre de décalage, et nouveauF est le score F du set de Ax 
				//choisi par la fonctionne phaseShift
			//shiftBy is the length of the shift (i.e. the number of positions that each
				//motif is shifted by), and nouveauF is the score F of the set of 
				//tDebuts (Ax) chosen by the phaseShift() function. 
			int shiftBy = (float)shiftNombres[0];
			float nouveauF = shiftNombres[1];

			//changer les debuts par le nombre shiftBy. 
			//change the start positions by the value of shiftBy
			for(int i = 0; i < nbSeq; i++) {
				unMotif.tDebut[i] += shiftBy;
			}

			//mettre à jour le tableau resultats avec les tDebuts changés. 
			//update the "resultats" table with the new tDebuts
			resultats[iteration][0] = iteration;
			for(int j = 0; j < nbSeq; j++) { 
				resultats[iteration][j+1] = unMotif.tDebut[j];
			}

			//guarder le vrai meilleur F. 
			//store the real best F score. 
			if(nouveauF > meilleur_F) {
				meilleur_F = nouveauF;
				index_de_meilleur_F = iteration;
			}
			
			//Si on ne change pas meilleur_index ni meilleur_F_actuel, l'algorithme reviendra au même maximum locale
				//qu'avant, donc il faux assurer que le nouveau meilleur_F_actuel soit le nouveauF trouvé par phaseShift. 
			//If neither meilleur_index nor meilleur_F_actuel change, the algorithm would return to the same local maximum
				//as before, so we need to ensure that the new meilleur_F_actuel is the nouveauF score found by the
				//phaseShift() function. 
			meilleur_index = iteration;
			meilleur_F_actuel = nouveauF;

			//on veut aussi plotter F et meilleur_F_actuel pour voir le performance de l'algorithme. On guarde les valeurs
				//de F et meilleur_F_actuel dans le table suivant:
			//We also want to plot F and meilleur_F_actuel to survey the performance of the algorithm. We keep the values 
				//of F and meilleur_F_actuel in the following table:
			suivre_F[iteration][0] = iteration;
			suivre_F[iteration][1] = nouveauF;
			suivre_F[iteration][2] = meilleur_F;

		} else {

			/*---------------------------------------------------------------------*/
			/*						OPTION 2: GIBBS SAMPLING					   */
			/*---------------------------------------------------------------------*/
			/*Dans cette option, on emploi le processus de Gibbs Sampling d'après Lawrence et al. 1993.*/ 
			/*In this option, we use the Gibbs Sampling process outlined in Lawrence et al. 1993.*/
			
			/*________________________________________________________________________________________*/
			/*PART 1: FIND A NEW ALIGNMENT BY ADJUSTING THE POSITION OF THE MOTIF IN A SINGLE SEQUENCE*/

			//définir un séquence aléatoirement à exlure.
			//Exclude a random sequence by chosing a number between 0 and the number of sequences randomly.
				//Note that the value of SEQ_A_EXCLURE is 1 + the index of the exluded sequence in lesSeq. 
				//The functions gibbs_functions.h account for this so that SEQ_A_EXCLURE can be inputted 
				//as an argument, but watch out! 
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

			float* rho = remplir_rho(unMotif, lesSeq, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											q(i,j)												 */
			/*-----------------------------------------------------------------------------------------------*/

			float** q = remplir_q(unMotif, C, rho, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											p(i)    											 */
			/*-----------------------------------------------------------------------------------------------*/

			float* p = remplir_p(unMotif, lesSeq, b, rho, nbSeq, SEQ_A_EXCLURE);

			/*-----------------------------------------------------------------------------------------------*/
			/*											Q(x)    											 */
			/*-----------------------------------------------------------------------------------------------*/

			//nMotif réprésent le nombre de motifs possible de longeur unMotif.lgMotif qui sont dans
				//la séquence exclue. 
			//nMotif represents the number of possible motifs of lengths unMotif.lgMotif which are in the 
				//excluded sequence.
			int nMotif = lesSeq[SEQ_A_EXCLURE-1].lg - unMotif.lgMotif + 1;

			float* Q = remplir_Q(nMotif, unMotif, lesSeq, q, SEQ_A_EXCLURE);

			/*-----------------------------------------------------------------------------------------------*/
			/*											P(x)    											 */
			/*-----------------------------------------------------------------------------------------------*/

			float* P = remplir_P(nMotif, unMotif, lesSeq, p, SEQ_A_EXCLURE);

			//Maintenant on calcule les ratios entre Q(x) et P(x) pour créer A(x). 
			//Now we calculate the ratios between Q(x) and P(x) to create A(x). 

			/*-----------------------------------------------------------------------------------------------*/
			/*											A(x)    											 */
			/*-----------------------------------------------------------------------------------------------*/
			
			float* A = remplir_A(Q, P, nMotif);

			
			//générer une numéro aléatoire pour choisir une des indices de A(x).
			//Generate a random number to choose one of the indices of A(x).
			double aleatoire = rand()/(double)(RAND_MAX);

			//guarder l'indece actuelle du motif dans la séquence exclue, pour les cas où l'algorithme choisit 
				//de ne pas guarder la nouvelle alignement. Si la nouvelle alignement n'est pas acceptée, on 
				//remettre l'indeceAncienne pour la position du motif dans la séquence exclue. 
			//Keep the current index of the motif in the excluded sequence, for the case where the l'algorithm 
				//choses not to accept the new alignment. If the new alignmnet is not accepted, we reset the 
				//position of the motif back to this value. 
			int indeceAncienne = unMotif.tDebut[SEQ_A_EXCLURE-1];

			//choisir une des indeces de A(x) par utiliser le variable "aleatoire", selon la distribution de 
				//A(x). L'indece choisit va être la nouvelle position du motif dans la séquence exclue, si
				//l'algorithme accepte la nouvelle alignment dans Partie 3. 
			//Chose one of the indices of A(x) by using the variable "aleatoire", according to the distribution
				//of A(x). The index chosen will be the new position of the motif in the excluded sequence, 
				//if the algorithm acceps the new alignment in part 3. 
			//Chose one of the indices 
			for(int i = 0; i < nMotif; i++) {
				if(A[i] > aleatoire) {
					unMotif.tDebut[SEQ_A_EXCLURE-1] = i;
					break;
				}
			}

			/*________________________________________________________________________________________*/
			/*PART 2: CALCULATE THE F-SCORE OF THE COMPLETE ALIGNMENT, INCLUDING THE EXCLUDED SEQUENCE*/

			//D'abord, il faut récalculer C(i,j), q(i,j), et p(j) avec la nouvelle alignement. C(i,j) et q(i,j)
				//doivent être avec l'alignement complête, autrement dit, avec la séquence exclue. 
			//First, we have to recalculate C(i,j) and q(i,j) with the new alignment. C(i,j) and q(i,j)
				//have to be with the complete alignment, in other words, with the excluded sequence. 

			/*-----------------------------------------------------------------------------------------------*/
			/*											C(i,j)												 */
			/*-----------------------------------------------------------------------------------------------*/

			//Cette fois, SEQ_A_EXCLURE est simplement une séquence qui n'existe pas, qui permet C(i,j) d'être 
				//sur l'alignement complête.
			//This time, SEQ_A_EXCLURE is simply a sequence which doesn't exist, which permits C(i,j) to be 
				//on the complete alignment. 
			C = remplir_C(unMotif, nbSeq+1, lesSeq, nbSeq);

			/*-----------------------------------------------------------------------------------------------*/
			/*											q(i,j)												 */
			/*-----------------------------------------------------------------------------------------------*/
			q = remplir_q(unMotif, C, rho, nbSeq);

			
			/*-----------------------------------------------------------------------------------------------*/
			/*											F    												 */
			/*-----------------------------------------------------------------------------------------------*/

			//calculer F
			//Calculate the F-score of the current alignment.
			float F = calculer_F(unMotif, C, q, p);


			//guarder le vrai meilleur F. 
			//If the new F is greater than the best F yet found (meilleur_F), then set meilleur_F to be 
				//equal to F and set the index_de_meilleur_F in the results table "resultats" to be equal
				//to the current iteration. 
			if(F > meilleur_F) {
				meilleur_F = F;
				index_de_meilleur_F = iteration;
			}

			/*________________________________________________________________________*/
			/*PART 2: SELECT OR REJECT THE THE NEW ALIGNMENT BASED ON SEVERAL CRITERIA*/
			
			/*--------------------------------------------*/
			/*CRITERE DE METROPOLIS (METROPOLIS CRITERION)*/
			/*--------------------------------------------*/
			//Pour ajouter plus d'hasard à l'algorithme, on inclut une critère de Metropolis, où on peut déplacer
				//le "meilleur_F_actuel" à un moindre score avec un probabilité qui est proportionel à la différence entre
				//le score F actuel et le score meilleur_F_actuel. La critère de Metropolis donne une deuxième façon 
				//d'échapper des maximums locals.
			//To add more randomness to the algorithm, we include a Metropollis criterion, where the algorithm can
				//move the "meilleur_F_actuel" score to a less score with a probability proportional to the difference between
				//the F score of the current alignment and the value of meilleur_F_actuel. This chance of chosing an alignment
				//with a less score gives another opportunity to escape local maximums.

			//"probabilite" décroisse quand F << meilleur_F_actuel.
			//"probabilite" decreases when F << meilleur_F_actuel.
			float probabilite = F/meilleur_F_actuel;
			double metropolis_chiffre = rand()/(double)(RAND_MAX);
			
			// resultats[iteration][0] = iteration;
			// for(int j = 0; j < nbSeq; j++) { 
			// 	resultats[iteration][j+1] = unMotif.tDebut[j];
			// }

			//Avec la critère de Metropolis, le variable "meilleur_index" va rester la même sauf si le nouveau 
				//alignement a un mieux score que le meilleur alignement, ou si l'alignement actuel est choisi
				//aléatoirement, selon la probabilité calculé dans le variable "probabilite". 
			//With the metropolis criterion, the "meilleur_index" variable will stay the same unless the 
				//new alignment has a better score than the best alignment ever found, or if the current
				//alignment is "selected" randomly accordning to a probability calculated in the "probabilite"
				//variable. 
			if(F > meilleur_F_actuel || probabilite >= metropolis_chiffre) {
				meilleur_F_actuel = F;
				meilleur_index = iteration;
			} else {
				//Si on n'accepte pas la nouvelle alignment, il faut remettre la position ancienne du motif dans 
					//la séquence exclue. 
				//If we don't accept the new alignment, we need to replace the former position of the motif 
					//in the excluded sequence. 
				unMotif.tDebut[SEQ_A_EXCLURE-1] = indeceAncienne;
			}

			//Ici, on guarde l'alignement guardé dans le "meilleur_index" de resultats. Le meilleur_index peut 
				//être l'index de l'itération actuelle, ou l'index d'une itération passée qui a trouvé la 
				//meilleur alignement. 
			//Here, we keep the alignment stored in the index indicated by the variable "meilleur_index" in the 
				//results. meilleur_index stores the index of the highest-score alignment in the "resultats"
				//table, which could be from some former iteration or from the current iterations. 
			
			// for(int col = 0; col < nbSeq; col++) {
			// 	unMotif.tDebut[col] = resultats[meilleur_index][col+1];
			// }

			//remplir "resultats" avec les tDebuts actuels
			//Fills "resultats" with the current tDebuts, or positions of the motifs. 
			resultats[iteration][0] = iteration;
			for(int j = 0; j < nbSeq; j++) { 
				resultats[iteration][j+1] = unMotif.tDebut[j];
			}

			//on veut aussi plotter F et meilleur_F_actuel pour voir le performance de l'algorithme. On guarde les valeurs
				//de F et meilleur_F_actuel dans le table suivant:
			//We also want to plot F and meilleur_F_actuel to see the performance of the algorithm. We keep the values of 
				//F and meilleur_F_actuel in the the table suivre_F. 
			suivre_F[iteration][0] = iteration;
			suivre_F[iteration][1] = meilleur_F_actuel;
			suivre_F[iteration][2] = meilleur_F;

			
			// //Calculate F without excluding a sequence so that we can track its progress:
			// calculer_F_pour_phaseShift(unMotif, lesSeq, numSeq, i)

			//vider la mémoire
			//Clear the memory
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
	//We kept the index of the best F score in "resultats". Now, the tDebut of the chosen motif
		//in each sequence is equal to the tDebuts in stored at at that index. 
	for(int col = 0; col < nbSeq; col++) {
		unMotif.tDebut[col] = resultats[index_de_meilleur_F][col+1];
	}
	
	/*-----------------------------------------------------------------------------------------------*/
	/*									AFFICHAGE (Printing)	         							 */
	/*-----------------------------------------------------------------------------------------------*/

	//afficher toutes les informations importants. 
	//Print all of the important information. 
	printf("Meilleur_F: %f \n", meilleur_F);
	printf("meilleur_F_actuel: %f \n", meilleur_F_actuel);
	printf("index_de_meilleur_F: %d \n", index_de_meilleur_F);
	for(int i=0; i < nbSeq; i++) {
		int tDebut = resultats[index_de_meilleur_F][i+1];
		char *seq = lesSeq[i].seq;
		printf("Sequence %d motif result: %d \n", i+1, tDebut);
		printf("%.*s\n", unMotif.lgMotif, seq + tDebut);
	}

	/*-----------------------------------------------------------------------------------------------*/
	/*											PLOTTING											 */
	/*-----------------------------------------------------------------------------------------------*/
	
	char * commandsForGnuplot[] = {"set title \"Alignments over all iterations\"", \
								   "plot 'data.temp' using 1:2 title 'Current F score' with lines, 'data.temp' using 1:3 title 'Best F score found' with lines"};
	for(int i=0;i<3;i++) {
		printf("%s\n", commandsForGnuplot[i]);
	}
    FILE * temp = fopen("data.temp", "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    for(int i = 0; i < NUM_ITERATIONS; i++) {
		fprintf(temp, "%-5f %-5f %-5f\n", suivre_F[i][0], suivre_F[i][1], suivre_F[i][2]);
		//fprintf(temp, "%-5f %-5f\n", suivre_F[i][0], suivre_F[i][2]);
	}

    for (i=0; i < NUM_COMMANDS; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }

	//imprimer le runtime
	//print the runtime
	end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Runtime: %f seconds\n", cpu_time_used);

	return 1;

}


