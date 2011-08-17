/*

Program implemented by Juan Manuel Lopez, Universidad Politecnica de Madrid
july-Augoust 2011
Implementation of Apodis 2

JET signal reading method provide by Jesus Vega (CIEMAT)


*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iso646.h>
#include "tipos.h"
#include "getfix/getfix.h"
#include "./tools/tools.h"
#include "./configuracion/config.h"



main(int argc, char *argv[]){
 
  char signalList[CHARMAX];  //input name of signals file
  int shotNumber;  // Input Shot number
  int Nsignals;    //Number the signal to use input+calculated

  signal *pSignal; //Pointer to struct signal, this struct contains all information of a signal

  FILE *fp;   //temporal pointer for file access
  int nlin;   //Counts the number of lines in file signallist. So the numebr of a raw signals
  char     txt[CHARMAX];   //Temporal array to hold oath names array for inout
  char   *ptr;   // temporal pointer to string
  char  title[53];
  char  units[11];


  int *t0;
  float *Maximums;
  float *Minimums;


  int i,j;  // for loops
  int error;  //Error indicator
  int ndata; //Number of datos of a signal


  /* Test of input parameters */
  if(argc != 3)
    {
      printf("\nProgram usage:\n\n");
      printf("> Main signalList shotNumber\n\n");
      printf("signalList: file with the JPF signal paths\n");
      printf("shotNumber: discharge to process\n");
      printf("\n\n");
      exit(0);
    }

  /* Input parameters to program variables */
  strcpy(signalList, argv[1]);
  shotNumber = atoi(argv[2]);


  read_config_files();

  printf ("Threshold %f \n",conf_Threshold);
  printf ("Sampling %f \n",conf_Sampling);
  printf ("Points  %d \n",conf_Npoints);
  printf ("Signals  %d \n",conf_Nsignals);
  printf ("Signal to calculate %d \n",conf_NcalculateSignals);
  printf ("NModels  %d \n",conf_NModels);
  printf ("Maximum %s \n",conf_PathMax);
  printf ("Minimum %s \n",conf_PathMin);
  printf ("Path d1 %s \n",conf_PathD1);
  printf ("path d2 %s \n",conf_PathD2);
  printf ("path d3 %s \n",conf_PathD3);
  printf ("path r %s \n",conf_PathR);

  Nsignals= conf_Nsignals + conf_NcalculateSignals; // Total number of signal to handle

  
  /* Number of signals to process */


  pSignal= (signal *) malloc(sizeof(signal)*Nsignals); //Alocate space for signals

  //Fill the struct with the path of the signal from the file signal list

  /* Number of signals to process */
  if((fp = fopen(signalList, "r")) == (FILE *) 0)
    {
      printf("\nError opening %s\n", signalList);
      free(pSignal);
      exit(0);
    }
  nlin = 0; /* Lines in the files (= number of signals) */
  while(fgets(txt, sizeof(txt), fp) != (char *) 0)
    {
      ptr = txt;
      while(*ptr > 32)
	++ptr;
      *ptr = 0;
      printf("Valor de txt %s \n",txt);

      strcpy((pSignal+nlin)->name,txt);
       ++nlin;
       if (nlin>conf_Nsignals){
	 printf("ERROR review Number of signal in Apodis2.cfg \n");
	 free(pSignal);
	 fclose(fp);
	 exit(0);
       }
    }
    fclose(fp);

    if (nlin != conf_Nsignals){
      	 printf("ERROR review Number of signal in Apodis2.cfg \n");
	 free(pSignal);
	 exit(0);
    }


  for(i=0;i<conf_Nsignals;i++){
    printf("Leido del struct %s \n",(pSignal+i)->name);
  }

  for (i=0;i<conf_Nsignals;i++){
    printf("Reading %s from shot %d\n",(pSignal+i)->name, shotNumber);
    fflush(stdout);

    getnwds_((pSignal+i)->name, &shotNumber, &ndata, &error, (long) strlen((pSignal+i)->name)); //Determine the space needed

    printf ("Debugg nData: %d error: %d\n", ndata, error);
    if ((!error) && (ndata > 0)){  //data are available
      if( ((pSignal+i)->pData = (float *) malloc(ndata*sizeof(float))) == (float *) 0 ){ //Allocate memory
              free(pSignal);
	      printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
	      exit(0);
       }
      if(((pSignal+i)->pTime = (float *) malloc(ndata*sizeof(float))) == (float *) 0){ //Allocate memory
	      free((pSignal+i)->pData);
              free(pSignal);
	      printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
	      exit(0);
       }

       // Pre set last char of returned text to 0 for c string  
       title [sizeof(title) - 1] = 0;
       units [sizeof(units) - 1] = 0;
       getdat_((pSignal+i)->name, &shotNumber, (pSignal+i)->pData, (pSignal+i)->pTime,&ndata , title, units, &error,(long) strlen((pSignal+i)->name), sizeof(title) - 1l, sizeof(units) - 1l);
       (pSignal+i)->nSamples= ndata;


    }
    else{
     
      printf("** Probable error due to non existing signal\n** shot %d, signal %s\n",shotNumber, (pSignal+i)->name);
      fflush(stdout);
      for(j = i; j <= 0; j--){
	   free((pSignal+i)->pData);
	   free((pSignal+i)->pTime);
      }
      free(pSignal);
      exit(0);

    }
  }

  //At this point the signal have been read form the JET BBDD

  for (i=1000;i<1016;i++){
    printf ("Signal t: %f valor %f \n",*((pSignal+1)->pTime+i),*((pSignal+1)->pData+i));

  }

  Maximums= (float *) malloc(sizeof(float)*Nsignals); //Alocate space for Maximums and Minimums
  Minimums= (float *) malloc(sizeof(float)*Nsignals); //Alocate space for Maximums and Minimums

  // Reading Maximums and Minimums values

  if((ReadFloatTxt(conf_PathMax, Maximums) != Nsignals) or (ReadFloatTxt(conf_PathMin, Minimums) !=  Nsignals)){
    printf ("Maximum or Minimus Files does not match with signal numbers \n");
    exit (0);
  }
    
 for(i=0;i<Nsignals;i++){
   (pSignal+i)->Max= *(Maximums+i);
   (pSignal+i)->Min= *(Minimums+i);
    printf("Leido del struct %f \n",(pSignal+i)->Max);
  }

 free(Maximums); //free array values are copied in the struct
 free(Minimums);

  
  //Look for t0 for the signals, the reference is Ipla signal 0


  t0= (int *) malloc(sizeof(int)*Nsignals); //Alocate space for initial times in the signals

  *t0= IndexEvent((pSignal)->pData,(pSignal)->nSamples,conf_Threshold,0); //Look for the time when Ipla < Threshold 

  // t0 for all signals, Only for raw signals

  for (i=1;i < conf_Nsignals;i++){
    *(t0+i)=  IndexEvent( (pSignal+i)->pTime,(pSignal+i)->nSamples,*((pSignal)->pTime+(*t0)),1);
  }

  //make pointer array and allocate memory for memory windows in models

  for (i=0;i<Nsignals;i++){ //Allocate array's pointer
     (pSignal+i)->pM= (float *) malloc(sizeof(float)*conf_NModels); //Allocate space for arrays pointers  
  }

  for (i=0;i<Nsignals;i++){ //Allocate memory buffers
     (pSignal+i)->pM= (float *) malloc(sizeof(float)*conf_NModels); //Allocate space 
     for(j=0;j<conf_NModels;j++){
       Maximums= (float *) malloc(conf_Npoints*sizeof(float)); //Allocate memory 
       
       *(((pSignal+i)->pM)+j)= &Maximums;

       
       //      if( ((pSignal+i)->(pM+j) = (float *) malloc(conf_Npoints*sizeof(float))) == (float *) 0 ){ //Allocate memory
       //       free(pSignal);
       //      printf("\n**** Not available memory for Model windows  %s, \n",(pSignal+i)->name);
       //	      exit(0);
       //}
     }
  }


  pSignal->Max= 45;
  (pSignal+1)->Max= 69;




}

