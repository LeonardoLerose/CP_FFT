/*
	This program computes the FFT in parallel using mpi
	Runs the program 'NUM_OF_EXEC' amount of times and computes the time
	for each iteration. Then calculate the average of the execution time
	Each process get a part (subtable) of the input complex vector (table).
	Then calculate even and odd for each subtable then P0 will gather these
	values and print them to the output file
	
COMPILE: mpcc FFT_Parallel.c -lm -o FFT_Parallel
RUN: llsubmit FFT_Parallel.job

Lerose Leonardo, Moretto Alberto, Perali Luca 
Universita degli studi di Padova
Laboratorio di Calcolo Parallelo
Implementazione in MPI della FFT
*/
#include <stdio.h>
#include <mpi.h>
#include <complex.h>
#include <math.h> /* for cos() and sin()*/
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265
#define INPUT_SIZE 16 /* Size of the input vector to create (e.g. 16384) */
#define ROW_LENGTH 50

int main()
{
	int my_rank, comm_sz;
	MPI_Init(NULL, NULL);					 /* start MPI */
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); /* how many processes are we using? */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /* which process is this? */
	double start, finish;
	double table[INPUT_SIZE][3];
	char row[ROW_LENGTH];
	char delim[] = " \n"; /*uso i delimitatori spazio e a capo per parsare il file di input*/
	FILE *outfile = NULL;
	FILE *inputfile = NULL;
	int h;


	if (my_rank == 0) /* if process 0, open outfile */
	{
		/*aggiunta creazione table da input
		TODO: fopen() deve prendere un input da argument 
		TODO: table non deve avere misura decisa a priori, va necessariamente creata dinamicamente*/
		inputfile = fopen("test1.txt", "rw+");
		if(inputfile == NULL)
		{
			printf("\n ERROR OPENING FILE, exit\n");
			return 1;
		}
		int i = 0;
		outfile = fopen("ParallelVersionOutput.txt", "w"); /* open from current directory */
		while (fgets(row, ROW_LENGTH, inputfile) != NULL)
		{
			char *ptr = strtok(row, delim);
			double index = atof(ptr);
			ptr = strtok(NULL, delim);
			double realPart = atof(ptr);
			ptr = strtok(NULL, delim);
			double imagPart = atof(ptr);
			table[i][0] = index;
			table[i][1] = realPart;
			table[i][2] = imagPart;
			/*printf("%.4f %.4f %.4f\n", table[i][0], table[i][1], table[i][2]);*/
			i++;
					
		}
		fclose(inputfile);
	}
	
	int i, k, n, j; /* Basic loop variables */

	double complex evenpart[(INPUT_SIZE / comm_sz / 2)];				 /*array to save the data for EVENHALF */
	double complex oddpart[(INPUT_SIZE / comm_sz / 2)];					 /*array to save the data for ODDHALF */
	double complex evenpartmaster[(INPUT_SIZE / comm_sz / 2) * comm_sz]; /* array to save the data for EVENHALF */
	double complex oddpartmaster[(INPUT_SIZE / comm_sz / 2) * comm_sz];  /* array to save the data for ODDHALF */
	double storeKsumreal[INPUT_SIZE];									 /* store the K real variable so we can abuse symmerty */
	double storeKsumimag[INPUT_SIZE];									 /* store the K imaginary variable so we can abuse symmerty */

	double subtable[(INPUT_SIZE / comm_sz)][3]; /* Each process owns a subtable from the table below  */

	int subtable_size = (INPUT_SIZE / comm_sz) * 3;														   /* how much to send and recieve */
	MPI_Scatter(table, subtable_size, MPI_DOUBLE, subtable, subtable_size, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* scatter the table to subtables */
	for (int i = 0; i<INPUT_SIZE / comm_sz; i++)
	{
		printf("procID: %d , subtable: %f \n", my_rank, subtable[i][0]);
	}
	for (k = 0; k < INPUT_SIZE / 2; k++) /* K coeffiencet Loop */
	{
		/* Variables used for the computation */
		double sumrealeven = 0.0; /* sum of real numbers for even */
		double sumimageven = 0.0; /* sum of imaginary numbers for even */
		double sumrealodd = 0.0;  /* sum of real numbers for odd */
		double sumimagodd = 0.0;  /* sum of imaginary numbers for odd */

		for (i = 0; i < (INPUT_SIZE / comm_sz) / 2; i++) /* Sigma loop EVEN and ODD */
		{
			double factoreven, factorodd = 0.0;
			int shiftevenonnonzeroP = subtable[2 * i][0];	/* used to shift index numbers for correct results for EVEN. */
			int shiftoddonnonzeroP =  subtable[2 * i + 1][0]; /* used to shift index numbers for correct results for ODD. */

			/* -------- EVEN PART -------- */
			double realeven = subtable[2 * i][1];						 /* Access table for real number at spot 2i */
			double complex imaginaryeven = subtable[2 * i][2];			 /* Access table for imaginary number at spot 2i */
			double complex componeeven = (realeven + imaginaryeven * I); /* Create the first component from table */
			if (my_rank == 0)											 /* if proc 0, dont use shiftevenonnonzeroP */
			{
				factoreven = ((2 * PI) * ((2 * i) * k)) / INPUT_SIZE; /*Calculates the even factor for Cos() and Sin() */
																	  /*   *********Reduces computational time********* */
			}
			else /*use shiftevenonnonzeroP */
			{
				factoreven = ((2 * PI) * ((shiftevenonnonzeroP)*k)) / INPUT_SIZE; /*Calculates the even factor for Cos() and Sin() */
																				  /*   *********Reduces computational time********* */
			}
			double complex comptwoeven = (cos(factoreven) - (sin(factoreven) * I)); /*Create the second component */

			evenpart[i] = (componeeven * comptwoeven); /*store in the evenpart array */

			/* -------- ODD PART -------- */
			double realodd = subtable[2 * i + 1][1];				  /*Access table for real number at spot 2i+1 */
			double complex imaginaryodd = subtable[2 * i + 1][2];	 /*Access table for imaginary number at spot 2i+1 */
			double complex componeodd = (realodd + imaginaryodd * I); /*Create the first component from table */
			if (my_rank == 0)										  /*if proc 0, dont use shiftoddonnonzeroP */
			{
				factorodd = ((2 * PI) * ((2 * i + 1) * k)) / INPUT_SIZE; /*Calculates the odd factor for Cos() and Sin()	*/
																		 /* *********Reduces computational time*********	 */
			}
			else /*use shiftoddonnonzeroP */
			{
				factorodd = ((2 * PI) * ((shiftoddonnonzeroP)*k)) / INPUT_SIZE; /*  Calculates the odd factor for Cos() and Sin() */
																				/* *********Reduces computational time********* */
			}

			double complex comptwoodd = (cos(factorodd) - (sin(factorodd) * I)); /*Create the second component */

			oddpart[i] = (componeodd * comptwoodd); /*store in the oddpart array */
		}

		/*Process ZERO gathers the even and odd part arrays and creates a evenpartmaster and oddpartmaster array*/
		MPI_Gather(evenpart, (INPUT_SIZE / comm_sz / 2), MPI_DOUBLE_COMPLEX, evenpartmaster, (INPUT_SIZE / comm_sz / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(oddpart, (INPUT_SIZE / comm_sz / 2), MPI_DOUBLE_COMPLEX, oddpartmaster, (INPUT_SIZE / comm_sz / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		if (my_rank == 0)
		{
			for (i = 0; i < (INPUT_SIZE / comm_sz / 2) * comm_sz; i++) /*loop to sum the EVEN and ODD parts */
			{
				sumrealeven += creal(evenpartmaster[i]); /*sums the realpart of the even half */
				sumimageven += cimag(evenpartmaster[i]); /*sums the imaginarypart of the even half */
				sumrealodd += creal(oddpartmaster[i]);   /*sums the realpart of the odd half */
				sumimagodd += cimag(oddpartmaster[i]);   /*sums the imaginary part of the odd half */
			}
			storeKsumreal[k] = sumrealeven + sumrealodd;				  /*add the calculated reals from even and odd */
			storeKsumimag[k] = sumimageven + sumimagodd;				  /*add the calculated imaginary from even and odd */
			storeKsumreal[k + INPUT_SIZE / 2] = sumrealeven - sumrealodd; /*ABUSE symmetry Xkreal + N/2 = Evenk - OddK */
			storeKsumimag[k + INPUT_SIZE / 2] = sumimageven - sumimagodd; /*ABUSE symmetry Xkimag + N/2 = Evenk - OddK */
			printf("%f \n", storeKsumreal[k]);
			if (k == 0)
			{
				fprintf(outfile, "\nTOTAL PROCESSED SAMPLES : %d\n", INPUT_SIZE);
			}
			fprintf(outfile, "FFT[%d]: %.4f + %.4fi \n", k, storeKsumreal[k], storeKsumimag[k]);

		}
	}

	if (my_rank == 0)
	{
		for(int k = INPUT_SIZE/2; k < INPUT_SIZE; k++)
		{
			fprintf(outfile, "FFT[%d]: %.4f + %.4f \n", k, storeKsumreal[k], storeKsumimag[k]);
		}
		fclose(outfile); /*CLOSE file ONLY proc 0 can. */
	}

	MPI_Barrier(MPI_COMM_WORLD); /*wait to all proccesses to catch up before finalize */
	MPI_Finalize();				 /*End MPI */
	return 0;
}
