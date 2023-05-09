// Ising2D T3.cpp
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#define N 13 // Trebuie pus cu 1 mai mult
#define JJ 1
#define T 3
#define Kb 1
#define NrPasiMonterCarlo 10001 // Trebuie pus cu 1 mai mult
#define _CRT_SECURE_NO_WARNINGS

int numaratorPozitivi = 0;
int numaratorNegativ = 0;
int spin[N][N];

int di[8] = {-1, 0, 1, 0, -1, -1, 1, 1};
int dj[8] = {0, 1, 0, -1, -1, 1, 1, -1};

// Define a recursive function to search for clusters
int search_cluster(int i, int j, int matrice[N][N], int vizitat[N][N])
{
	// Mark the site as visited
	vizitat[i][j] = 1;
	// Check the neighboring sites

	for (int k = 0; k < 8; k++)
	{
		int ni = i + di[k];
		int nj = j + dj[k];
		if (matrice[ni][nj] != 0 && vizitat[ni][nj] == 0 && matrice[ni][nj] == matrice[i][j])
		{
			// Recursively search the neighboring site if it has the same spin and has not been visited
			return 1 + search_cluster(ni, nj, matrice, vizitat);
		}
	}

	return 1;
}

void medii_clustere(int a[N][N], FILE *fptr)
{

	int viz[N][N] = {0};

	int nr_clustere_poz = 0;
	int nr_clustere_neg = 0;
	int sum_marimi_poz = 0;
	int sum_marimi_neg = 0;

	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			if (!viz[i][j])
			{
				int marime_cluster = search_cluster(i, j, a, viz);
				if (a[i][j] == 1)
				{
					nr_clustere_poz++;
					sum_marimi_poz += marime_cluster;
				}
				else
				{
					nr_clustere_neg++;
					sum_marimi_neg += marime_cluster;
				}
			}
		}
	}

	float average_poz = (float)sum_marimi_poz / (float)nr_clustere_poz;
	float average_neg = (float)sum_marimi_neg / (float)nr_clustere_neg;
	fprintf(fptr, "Nr. clustere poz: %d \t Medie clustere pozitive: %f \n", nr_clustere_poz, average_poz);
	fprintf(fptr, "Nr. clustere neg: %d \t Medie clustere negative: %f \n", nr_clustere_neg, average_neg);
}

/*int corelatie(int spin[N][N])
{
	FILE *fptr;
	fptr = fopen("C.txt", "w");
	int C = 0;
	int Suma = 0;
	for (int x = 1; x < N - 1; x++)
	{
		for (int y = 1; y < N - 1; y++)
		{
			Suma = spin[x - 1][y] + spin[x + 1][y] + spin[x][y - 1] + spin[x][y + 1];
			C += Suma * spin[x][y];
			// Suma = 0;
		}
	}
	printf("%d", C);
	fprintf(fptr, "%d", C);
	fclose(fptr);
	return C;
}*/

void Metropolis()
{
	// Matricea
	int i, j;
	int Sv = 0;						// Suma spinilor vecini - Nord Sud Est Vest
	int M[NrPasiMonterCarlo] = {0}; // Magnetizatia
	double P;
	double r;
	int dE = 0;
	FILE *fptr;

	for (i = 0; i < N; i++) // Atribuim tuturor spinilor 1
	{
		for (j = 0; j < N; j++)
		{
			spin[i][j] = 1;
		}
	}

	for (j = 0; j < N; j++) // Egalam marginile cu zero
	{
		spin[0][j] = 0;
		spin[N - 1][j] = 0;
		spin[j][0] = 0;
		spin[j][N - 1] = 0;
	}

	int k = 0;
	srand(time(NULL));
	FILE *rezultat_ptr = fopen("Rezultat.txt", "w");
	FILE *imagini_ptr = fopen("ImaginiIntermediare.txt", "w");
	for (int pasiMC = 0; pasiMC < NrPasiMonterCarlo; pasiMC++) // Pasi Monte Carlo
	{
		for (int contor1 = 0; contor1 < N * N; contor1++)
		{
			int x = int(rand() / (RAND_MAX + 1.0) * N); // Alegem spinii carora vrem sa le schimbam semnul
			int y = int(rand() / (RAND_MAX + 1.0) * N);

			if (x >= 1 && x <= N - 2 && y >= 1 && y <= N - 2) // Fara margini
			{
				Sv = spin[x - 1][y] + spin[x + 1][y] + spin[x][y - 1] + spin[x][y + 1]; // Calculam diferenta de energie
				dE = 2 * JJ * Sv * spin[x][y];
				P = exp(-dE / T); // Probabilitatea ca spinii sa se schimbe	Beta =
				r = (rand() / (RAND_MAX + 1.0));
				if (r < P)
				{
					spin[x][y] = -spin[x][y];
				}
			}
		}
		fprintf(rezultat_ptr, "\nIteratia numarul:%d\n", pasiMC);
		medii_clustere(spin, rezultat_ptr);
		// sa printeze din 2000 in 2000
		if (pasiMC % 2000 == 0)
		{
			fprintf(imagini_ptr, "Iteratia nr: %d\n", pasiMC);
			for (int linie = 0; linie < N; linie++)
			{
				for (int coloana = 0; coloana < N; coloana++)
				{

					fprintf(imagini_ptr, "%2d\t%2d\t%2d\n", linie, coloana, spin[i][j]);
				}
			}
			fprintf(imagini_ptr, "\n");
		}
	}
	fclose(rezultat_ptr);
	// scrie matricea spin in fisier
	fptr = fopen("Matrice.txt", "w");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{

			fprintf(fptr, "%2d\t%2d\t%2d\n ", i, j, spin[i][j]);
		}
	}
	fclose(fptr);
	// afiseara matricea spin pe ecran
	printf("\n");
	/*for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
		printf("\n");
	}*/
}

int main()
{

	Metropolis();
	for (int i = 0; i < N; i++)

		for (int j = 0; j < N; j++)
		{
			printf("%2d ", spin[i][j]);
		}
	printf("\n");

	// numara spinii de acelasi semn aflati unu langa celalalt pe aceeasi linie
	// numarare(spin);

	// calculeaza factorul de corelatie
	// corelatie(spin);

	return 0;
}