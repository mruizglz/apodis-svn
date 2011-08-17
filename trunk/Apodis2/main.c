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

void PRINTSIGNAL( signal *entrada);


main(int argc, char *argv[]){
 
  char signalList[CHARMAX];  //input name of signals file
  int shotNumber;  // Input Shot number
  int Nsignals;    //Number the signal to use input+calculated

  signal *pSignal; //Pointer to struct signal, this struct contains all information of a signal

  FILE *fp;   //temporal pointer for file access
  int nlin;   //Counts the number of lines in file signallist. So the numebr of a raw signals
  char     txt[CHARMAX];   //Temporal array to hold path names array for input
  char   *ptr;   // temporal pointer to string
  char  title[53];
  char  units[11];


  int *t0;
  float *Maximums;
  float *Minimums;


  int i,j,t,z;  // for loops
  int error;  //Error indicator
  int ndata; //Number of datos of a signal

  float **pArray; //pointer to array
  float *pBuffer; //Pointer too buffer
  float *pBuffer1; //Pointer too buffer
  float *pBuffer2; //Pointer too buffer

  int  *ToNormalize; // pointer to array with actions about normalize proccess
                     //[12]={TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE};
  char **PathMx;  //Path to files M for model

  model *pModel;  //Pointer to struct model

  float dummy;
  float **Mdummy;

  int **rowptr;

  float R[4];   //Este hay que pasarlo tambien a dinamico
  float ModelParts[11][3];//Este hay que pasarlo tambien a dinamico

  double (*infft)[2];   //Buffer for FFT 
  double (*outfft)[2];

  double resabs[32];



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
  printf ("path Model %s \n",conf_PathModel);
  printf ("path r %s \n",conf_PathR);
  printf ("path Normalize  %s \n",conf_PathNormalize);

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
       (pSignal+i)->Npoints= conf_Npoints;

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
  /*
  for (i=1143;i<1155;i++){
    printf ("Signal t: %f valor %f \n",*((pSignal+0)->pTime+i),*((pSignal+0)->pData+i));

  }
  */


  Maximums= (float *) malloc(sizeof(float)*Nsignals); //Alocate space for Maximums and Minimums
  Minimums= (float *) malloc(sizeof(float)*Nsignals); //Alocate space for Maximums and Minimums
  ToNormalize= (int *) malloc(sizeof(int)*Nsignals); //Alocate space for normaliza Flag

  // Reading Maximums and Minimums values

  if((ReadFloatTxt(conf_PathMax, Maximums) != Nsignals) or (ReadFloatTxt(conf_PathMin, Minimums) !=  Nsignals or (ReadNormalizeTxt(conf_PathNormalize, ToNormalize) != Nsignals))){
    printf ("Maximum, Minimus, or Normalize  Files does not match with signal numbers \n");
    exit (0);
  }
    
 for(i=0;i<Nsignals;i++){
   (pSignal+i)->Max= *(Maximums+i);
   (pSignal+i)->Min= *(Minimums+i);
   (pSignal+i)->Normalize= *(ToNormalize+i);   
   //   printf("Leido del struct %d \n",(pSignal+i)->Normalize);
  }

 free(Maximums); //free array values are copied in the struct
 free(Minimums);
 free(ToNormalize);

  
  //Look for t0 for the signals, the reference is Ipla signal 0


  t0= (int *) malloc(sizeof(int)*Nsignals); //Alocate space for initial times in the signals

  *t0= IndexEvent((pSignal)->pData,(pSignal)->nSamples,conf_Threshold,0); //Look for the time when Ipla < Threshold 

  printf("Indice de t0 %d \n",*t0);

  // t0 for all signals, Only for raw signals

  for (i=1;i < conf_Nsignals;i++){
    *(t0+i)=  IndexEvent( (pSignal+i)->pTime,(pSignal+i)->nSamples,*((pSignal)->pTime+(*t0)),1);
  }
  
  //make pointer array and allocate memory for memory windows in models

  for (i=0;i<Nsignals;i++){ //Allocate array's pointer and pointer to resamplig time
    (pSignal+i)->pM=  malloc(sizeof(float)*conf_NModels); //Allocate space for arrays pointers  
  }
 
  for (i=0;i<Nsignals;i++){ //Allocate memory buffers
    (pSignal+i)->pTimeR= malloc(sizeof(float)*conf_Npoints); //Buffer for resampling time
    for(j=0;j<conf_NModels;j++){
      *((pSignal+i)->pM+j)= (float *) malloc(conf_Npoints*sizeof(float));
      if(*((pSignal+i)->pM+j) == NULL ){ //Do so for easy code ready
	free(pSignal);
        printf("\n**** Not available memory for Model windows  %s, \n",(pSignal+i)->name );
        exit(0);
      }
    }
  }
 
  //Here the raw signal are read, and all buffers are allocated for all signal
  //Now the buffer are going to be filled with resampling data window.

  
  for (i=0;i<conf_Nsignals;i++){ //Resampling for rae signals only 
    //     printf ("struct status signal %s : \n",(pSignal+i)->name);
    // printf ("maximum: %f minimum: %f Npuntos: %d nsamples: %d \n",(pSignal+i)->Max,(pSignal+i)->Min,(pSignal+i)->Npoints,(pSignal+i)->nSamples);
     // printf ("pData: %d pTime: %d pTimeR: %d \n",(pSignal+i)->pData,(pSignal+i)->pTime,(pSignal+i)->pTimeR);
     for (j=conf_NModels-1;j>=0;j--){
       // printf ("pM %d : ",(pSignal+i)->pM+j);
        t=resampling ( *(t0+i), conf_Sampling,*((pSignal+i)->pM+j) ,(pSignal+i));
       *(t0+i)+=t;
       printf("\n Cambio senal \n");
       pBuffer=*((pSignal+i)->pM+j);
       //       for(z=0;z<32;z++)
       //	 printf("Valor normalizado: %f ",*(pBuffer+z));
     }
     //   printf("\n");


  }


  //Here the raw signal are reading and ready to be used
 
  //Now start the caculate signals. Signal number 8 (from number 3) the difference signal
  //Signal 9 (from 6 and 7). Signal 6 is adjust to a minimum value of 1000 and then is divided 
  //for signal 7.

  //Siganl 8. Index=7 (8 - 1)
  //Make the new signal x[n]=x[n+1]-x[n], the last point is the previous point

  for (j=conf_NModels-1;j>0;j--){  //two parts loop for use the last loop M3 to generate the timeR
     // printf ("pM %d : ",(pSignal+i)->pM+j);
    pBuffer= *((pSignal+7)->pM+j);  //Used for easy code reading
    pBuffer1= *((pSignal+2)->pM+j);
    for (i=0;i<conf_Npoints-1;i++){
      *(pBuffer+i)= normalize ((pSignal+7)->Max, (pSignal+7)->Min, *(pBuffer1+i+1)-*(pBuffer1+i));
    }
    *(pBuffer+i)= *(pBuffer+i-1);
  }
  
  for (j=0;j>=0;j--){                           //genarate M3 timeR and data
     // printf ("pM %d : ",(pSignal+i)->pM+j);
    pBuffer= *((pSignal+7)->pM+j);  //Used for easy code reading
    pBuffer1= *((pSignal+2)->pM+j);
   
    for (i=0;i<conf_Npoints-1;i++){
      *(pBuffer+i)= normalize ((pSignal+7)->Max, (pSignal+7)->Min, *(pBuffer1+i+1)-*(pBuffer1+i));
      *((pSignal+7)->pTimeR+i)= *((pSignal+2)->pTimeR+i);
    }
    *(pBuffer+i)= *(pBuffer+i-1);
    *((pSignal+7)->pTimeR+i)= *((pSignal+7)->pTimeR+i-1);
  }
        
  //Signal 8. Done

  //Signal 9. difference between signals 6 and 7. The signal is number 8 (9-1)

  for (j=conf_NModels-1;j>0;j--){  //two parts loop for use the last loop M3 to generate the timeR
     // printf ("pM %d : ",(pSignal+i)->pM+j);
    pBuffer= *((pSignal+8)->pM+j);  //Used for easy code reading
    pBuffer1= *((pSignal+5)->pM+j);
    pBuffer2= *((pSignal+6)->pM+j);
    for (i=0;i<conf_Npoints;i++){
      *(pBuffer+i)= normalize ((pSignal+8)->Max, (pSignal+8)->Min, (*(pBuffer1+i) < 1000 ? 1000 :*(pBuffer1+i)) / *(pBuffer2+i)); //Ojo 1000 a configuracion
    }
  }

  
  for (j=0;j>=0;j--){                           //genarate M3 timeR and data
     // printf ("pM %d : ",(pSignal+i)->pM+j);
    pBuffer= *((pSignal+8)->pM+j);  //Used for easy code reading
    pBuffer1= *((pSignal+5)->pM+j);
    pBuffer2= *((pSignal+6)->pM+j);
   
    for (i=0;i<conf_Npoints;i++){
      *(pBuffer+i)= normalize ((pSignal+8)->Max, (pSignal+8)->Min, (*(pBuffer1+i) < 1000 ? 1000 :*(pBuffer1+i)) / *(pBuffer2+i)); //Ojo 1000 a configuracion
      *((pSignal+8)->pTimeR+i)= *((pSignal+5)->pTimeR+i);
    }
  }



  // PRINTSIGNAL (pSignal);
  //  PRINTSIGNAL (pSignal+7);


  //Reading Model Information

   PathMx=  malloc(sizeof(char *)*conf_NModels); //Alocate space for Path to Models properties
 
   for (i=0;i<conf_NModels;i++){ //Allocate array's pointer 
     *(PathMx+i)= (char *) malloc(sizeof(char)*CHARMAX); //Allocate space for arrays
   }
   if (ReadModelTxt(conf_PathModel, PathMx) != conf_NModels){ //Check the number of models
     printf ("Review Lenght FILE Model description \n");
     exit(0);
   }

   pModel= (model *) malloc(sizeof(model)*conf_NModels); //Alocate space for signals
   if (pModel == NULL){
       free(pSignal);
       printf ("Not space for Mpdel struct \n");
       exit(0);
   }

   for (i=0;i<conf_NModels;i++){ //Allocate array data value for model 
     N_vectors(*(PathMx+i),&(pModel+i)->nvectors , &(pModel+i)->coef_vector);
     printf("Tamanos Modelo: %d vectores: %d coeficientes: %d \n",i,(pModel+i)->nvectors , (pModel+i)->coef_vector);

     (pModel+i)->alfa = malloc(sizeof(float )*(pModel+i)->nvectors); //Allocate space for arrays
     if (Mdummy == NULL){
        free(pSignal);
        printf ("Not space for alfa data array \n");
        exit(0);
     }

     Mdummy = malloc(sizeof(float *)*(pModel+i)->nvectors); //Allocate space for arrays of data

     if (Mdummy == NULL){
       free(pSignal);
       printf ("Not space for Mmodel data array \n");
       exit(0);
     }
     for (j=0;j<(pModel+i)->nvectors;j++){
       Mdummy[j]= malloc(sizeof(float)*(pModel+i)->coef_vector); //Allocate space for arrays
       if (Mdummy[j] == NULL){
	 free(pSignal);
	 printf ("Not space for Mmodel data array \n");
	 exit(0);
       }

     }

     (pModel+i)->data= Mdummy;

     M_values (*(PathMx+i), pModel+i);
   }

   Read_M (conf_PathR, R);

   //here the model is read

   //calculamos los valores del medelo de las tres ventanas
   //esto hay que meterlo dentro de la estructira de punteros


   for (i=0;i<conf_NModels;i++){ //mean ipla 
     ModelParts[0][i]= Mean (*((pSignal)->pM+i), conf_Npoints);
   }

   for (i=0;i<conf_NModels;i++){ //des ipla 
     

     fft(conf_Npoints, infft, outfft);  
 
     Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT
 
     ModelParts[0][i]= Desv (resabs, conf_Npoints/2);

   }


  //    N_vectors("p.txt" , &dmmy, &dmmy1);

  //   printf ("resultados %d %d \n",dmmy,dmmy1);


}





void PRINTSIGNAL( signal *entrada){
  int i;

  printf ("Valores comunes Max: %f  Min: %f Pts: %d \n",entrada->Max, entrada->Min, entrada->Npoints);
  printf ("ind \t Tiempo \t\t valor M3 \t valor M2 \t valor M1 \n");
  for (i=0;i<entrada->Npoints;i++){
    printf ("%d \t tr: %f  \t  %f  \t  %f \t  %f \n", i, *(entrada->pTimeR+i),*(entrada->pM+2+i),*(entrada->pM+1+i),*(entrada->pM+i));
  }

}
