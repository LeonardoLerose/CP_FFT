/*
Lerose Leonardo, Moretto Alberto, Perali Luca 
Università degli studi di Padova
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
#define ROW_LENGTH 50 /* Used in parsing inputFile */

int main(int argc, char* argv[])
{		
	if(argc != 3)
	{
		printf(" ERROR, BAD ARGUMENTS");
		return 1;
	}

	int my_rank, comm_sz; /* Rank of proc, Number of proc */
	char row[ROW_LENGTH]; /* A single line of file */
	char delim[] = " \n"; /* " " and "\n" used to parse the input */
	
	FILE *inputfile = NULL; /* File that contains the input vector, choosen with argv[2] */
	FILE *outfile = NULL;	/* File that contains the FFT of input vector */ 


	char* vectorSize = argv[1];		/* Size of the vector of which perform FFT */
	char* inputFileName = argv[2];	/* Name of the file where is stored the input vector */
	
	int dimensions = atoi(vectorSize);
	double table[dimensions][3];	/* Create a table with <dimensions> rows and 3 columns: column0: index, column1: Re, column2: Im */

	MPI_Init(NULL, NULL);					 
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	if (my_rank == 0) 
	{
		inputfile = fopen(inputFileName, "r");
		if(inputfile == NULL)
		{
			printf("ERROR, OPENING FILE");
			return 1;
		}

		outfile = fopen("ParallelVersionOutput.txt", "w"); 	/* DA RIMUOVERE, Usiamo outfile quello specificato nel file .job (usato per comodità in questo momento) */

		int i = 0;
		while (fgets(row, ROW_LENGTH, inputfile) != NULL)	/* Parse inputfile one row at the time */
		{
			char *next = strtok(row, delim);
			double index = atof(next);
			next = strtok(NULL, delim);
			double realPart = atof(next);
			next = strtok(NULL, delim);
			double imaginaryPart = atof(next);
			table[i][0] = index;	/*populate the table */
			table[i][1] = realPart;
			table[i][2] = imaginaryPart;
			i++;
					
		}
		fclose(inputfile);
	}
	
	int i, k, n, j;		/* Vedere se metterli dentro al for */

	double complex evenpart[(dimensions / comm_sz / 2)];				 /* used to store the even part of data in the single proc */ 	/* nome = partE */
	double complex oddpart[(dimensions / comm_sz / 2)];					 /* used to store the odd part of data in the single proc */	/* nome = partO */
	double complex evenpartmaster[(dimensions / comm_sz / 2) * comm_sz]; /* used to store the even part of data in the master  */		/* nome = partEM */
	double complex oddpartmaster[(dimensions / comm_sz / 2) * comm_sz];  /* used to store the odd part of data in the master */			/* nome = partOM */
	double storeKsumreal[dimensions];									 /* Here will be saved the real values of the output FFT vector */ /* nome = sumKRe */
	double storeKsumimag[dimensions];									 /* Here will be saved the imaginary of the output FFT vector */	/* nome = sumKIm */
	double subtable[(dimensions / comm_sz)][3];							 /* Part of table owned by a single process  */
	int subtable_size = (dimensions / comm_sz) * 3;						 /* How big each subtable */
	
	MPI_Scatter(table, subtable_size, MPI_DOUBLE, subtable, subtable_size, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* Send data in table to subtables */

	for (int i = 0; i<dimensions / comm_sz; i++)
	{
		printf("procID: %d , subtable: %f \n", my_rank, subtable[i][0]);
	}
	for (k = 0; k < dimensions / 2; k++) /* K coeffiencet Loop */
	{
		/* Variables used for the computation */
		double sumrealeven = 0.0; /* sum of real numbers for even */
		double sumimageven = 0.0; /* sum of imaginary numbers for even */
		double sumrealodd = 0.0;  /* sum of real numbers for odd */
		double sumimagodd = 0.0;  /* sum of imaginary numbers for odd */

		for (i = 0; i < (dimensions / comm_sz) / 2; i++) /* Sigma loop EVEN and ODD */
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
				factoreven = ((2 * PI) * ((2 * i) * k)) / dimensions; /*Calculates the even factor for Cos() and Sin() */
																	  /*   *********Reduces computational time********* */
			}
			else /*use shiftevenonnonzeroP */
			{
				factoreven = ((2 * PI) * ((shiftevenonnonzeroP)*k)) / dimensions; /*Calculates the even factor for Cos() and Sin() */
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
				factorodd = ((2 * PI) * ((2 * i + 1) * k)) / dimensions; /*Calculates the odd factor for Cos() and Sin()	*/
																		 /* *********Reduces computational time*********	 */
			}
			else /*use shiftoddonnonzeroP */
			{
				factorodd = ((2 * PI) * ((shiftoddonnonzeroP)*k)) / dimensions; /*  Calculates the odd factor for Cos() and Sin() */
																				/* *********Reduces computational time********* */
			}

			double complex comptwoodd = (cos(factorodd) - (sin(factorodd) * I)); /*Create the second component */

			oddpart[i] = (componeodd * comptwoodd); /*store in the oddpart array */
		}

		/*Process ZERO gathers the even and odd part arrays and creates a evenpartmaster and oddpartmaster array*/
		MPI_Gather(evenpart, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, evenpartmaster, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(oddpart, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, oddpartmaster, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		if (my_rank == 0)
		{
			for (i = 0; i < dimensions / 2; i++) /*loop to sum the EVEN and ODD parts */
			{
				sumrealeven += creal(evenpartmaster[i]); /*sums the realpart of the even half */
				sumimageven += cimag(evenpartmaster[i]); /*sums the imaginarypart of the even half */
				sumrealodd += creal(oddpartmaster[i]);   /*sums the realpart of the odd half */
				sumimagodd += cimag(oddpartmaster[i]);   /*sums the imaginary part of the odd half */
			}
			storeKsumreal[k] = sumrealeven + sumrealodd;				  /*add the calculated reals from even and odd */
			storeKsumimag[k] = sumimageven + sumimagodd;				  /*add the calculated imaginary from even and odd */
			storeKsumreal[k + dimensions / 2] = sumrealeven - sumrealodd; /*ABUSE symmetry Xkreal + N/2 = Evenk - OddK */
			storeKsumimag[k + dimensions / 2] = sumimageven - sumimagodd; /*ABUSE symmetry Xkimag + N/2 = Evenk - OddK */
			printf("%f \n", storeKsumreal[k]);
			if (k == 0)
			{
				fprintf(outfile, "\nTOTAL PROCESSED SAMPLES : %d\n", dimensions);
			}
			fprintf(outfile, "FFT[%d]: %.2f + %.2fi \n", k, storeKsumreal[k], storeKsumimag[k]);

		}
	}
	
	if (my_rank == 0)
	{
		for(int k = dimensions/2; k < dimensions; k++)
		{
			fprintf(outfile, "FFT[%d]: %.2f + %.2fi \n", k, storeKsumreal[k], storeKsumimag[k]);
		}
		fclose(outfile); 
	}

	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize();				 
	return 0;
}
