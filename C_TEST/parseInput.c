#include<stdio.h>
#include<string.h> 
#include <stdlib.h>
 
int main(void) 
{
    FILE* fd = NULL; 
    char ch;
    char row[50];
    char delim[] = " \n"; //list of delimit
    fd = fopen("test.txt","rw+"); 
    if(NULL == fd)
    { 
        printf("\n fopen() Error!!!\n"); 
        return 1;
    }
    int i = 0;
    //crea struttura per memorizzare dati
    double table[10][3];      
    //leggi i dati, uno per volta e inserisci nell'array
    while ((fgets(row, 50, fd)) != NULL)
    {
        //first number
        char* ptr = strtok(row, delim);
        // printf("%s\n", ptr);
        double n = atof(ptr);
        //second number
        ptr = strtok(NULL, delim);
        // printf("%s\n", ptr);
        double m = atof(ptr);
        table[i][1] = n;
        table[i][2] = m;
        i++;
    }
    for (int i = 0; i < 10; i++)
    {
        printf("Re: %f , Im: %f \n", table[i][1], table[i][2]);
    }
    
    
    fclose(fd);

    return 0; 
}


