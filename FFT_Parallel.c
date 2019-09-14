/*
Lerose Leonardo, Moretto Alberto, Perali Luca 
Università degli studi di Padova
Laboratorio di Calcolo Parallelo
Implementazione in MPI della FFT
*/

#include <stdio.h>
#include <mpi.h>
#include <complex.h>
#include <math.h> 
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265
#define ROW_LENGTH 50 /* Used in parsing inputFile */

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		perror(" ERROR, BAD ARGUMENTS");
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

	double tot_time = -MPI_Wtime();		
	double loadinput_time;
    double scatter_time;
    double evenodd_time;
    double gather_time;

	if (my_rank == 0) 
	{
        loadinput_time -= MPI_Wtime();
		inputfile = fopen(inputFileName, "r");
		if(inputfile == NULL)
		{
			perror("ERROR, OPENING FILE");
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
        loadinput_time += MPI_Wtime();
	}
	
	int i, k, n, j;		/* VEDERE SE METTERLI DENTRO AL FOR */

	double complex partEven[(dimensions / comm_sz / 2)];				 /* used to store the even part of data in the single proc */ 	 
	double complex partOdd[(dimensions / comm_sz / 2)];					 /* used to store the odd part of data in the single proc */	
	double complex partEvenmaster[(dimensions / comm_sz / 2) * comm_sz]; /* used to store the even part of data in the master  */		
	double complex partOddmaster[(dimensions / comm_sz / 2) * comm_sz];  /* used to store the odd part of data in the master */			
	double sumKRe[dimensions];									 		 /* Here will be saved the real values of the output FFT vector */
	double sumKIm[dimensions];									 		 /* Here will be saved the imaginary of the output FFT vector */
	double subtable[(dimensions / comm_sz)][3];							 /* Part of table handled by a single process  */
	int subtable_size = (dimensions / comm_sz) * 3;						 /* How big each subtable is */

    scatter_time -= MPI_Wtime();
	MPI_Scatter(table, subtable_size, MPI_DOUBLE, subtable, subtable_size, MPI_DOUBLE, 0, MPI_COMM_WORLD); /* Send data in table to subtables */
    scatter_time += MPI_Wtime();

    evenodd_time -= MPI_Wtime();
	for (k = 0; k < dimensions / 2; k++) /* 0 to input.size/2 loop */
	{
		/* Variables used for the computation */
		double sumReE = 0.0; /* sum of real numbers for even */
		double sumImE = 0.0; /* sum of imaginary numbers for even */
		double sumReO = 0.0;  /* sum of real numbers for odd */
		double sumImO = 0.0;  /* sum of imaginary numbers for odd */

		for (i = 0; i < (dimensions / comm_sz) / 2; i++) /* Sigma loop EVEN and ODD */
		{
			double factorEven, factorOdd = 0.0;
			int shiftE = subtable[2 * i][0];			/* shift the right amount proc that have rank !=0 , Even part */
			int shiftO =  subtable[2 * i + 1][0]; 	/* shift the right amount proc that have rank !=0 , Odd part */

			double reE = subtable[2 * i][1];					/* real number at position 2i */
			double complex imE = subtable[2 * i][2];			/* imaginary number at position 2i */
			double complex firstCompE = (reE + imE * I); 		
			
			double reO = subtable[2 * i + 1][1];				/* real number at position 2i+1 */
			double complex imO = subtable[2 * i + 1][2];	 	/* imaginary number at position 2i+1 */
			double complex firstCompO = (reO + imO * I); 

			factorEven = ((2 * PI) * ((shiftE)*k)) / dimensions; 
			factorOdd = ((2 * PI) * ((shiftO)*k)) / dimensions; 	

			double complex secondCompE = (cos(factorEven) - (sin(factorEven) * I)); /*Create the second component */
			double complex secondCompO = (cos(factorOdd) - (sin(factorOdd) * I)); 	

			partEven[i] = (firstCompE * secondCompE); 	/*store in the partEven array */
			partOdd[i] = (firstCompO * secondCompO); 	/*store in the partOdd array */

		}

		
        gather_time -= MPI_Wtime();
		/* populate the master array */
		MPI_Gather(partEven, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, partEvenmaster, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(partOdd, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, partOddmaster, (dimensions / comm_sz / 2), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
        gather_time += MPI_Wtime();

		if (my_rank == 0)
		{
			for (i = 0; i < dimensions / 2; i++) /*loop to sum the EVEN and ODD parts */
			{
				sumReE += creal(partEvenmaster[i]);	/*sum of Re even part */
				sumImE += cimag(partEvenmaster[i]);	/*sum of Im even part */
				sumReO += creal(partOddmaster[i]);	/*sum of Re odd part */
				sumImO += cimag(partOddmaster[i]);	/*sum of Im odd part */
			}

			sumKRe[k] = sumReE + sumReO;				  
			sumKIm[k] = sumImE + sumImO;				  
			sumKRe[k + dimensions / 2] = sumReE - sumReO; 
			sumKIm[k + dimensions / 2] = sumImE - sumImO;
			 
			if (k == 0)
			{
				fprintf(outfile, "TOTAL PROCESSED SAMPLES : %d\n", dimensions);
			}
			fprintf(outfile, "FFT[%d]: %.2f + %.2fi \n", k, sumKRe[k], sumKIm[k]);

		}
	} /* 0 to input.size/2 loop */
    evenodd_time += MPI_Wtime();

	if (my_rank == 0)
	{
		for(int k = dimensions/2; k < dimensions; k++)
		{
			fprintf(outfile, "FFT[%d]: %.2f + %.2fi \n", k, sumKRe[k], sumKIm[k]);
		}
		fclose(outfile); 
	}


    tot_time += MPI_Wtime();
    printf("%d,%d,%f,%f,%f,%f\n", my_rank, comm_sz, tot_time, scatter_time, evenodd_time, gather_time);
    fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD); /*wait to all proccesses to catch up before finalize */
	MPI_Finalize();				 /*End MPI */
	return 0;
}
