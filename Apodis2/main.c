
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
#include "./fft/fft.h"

void PRINTSIGNAL( signal *entrada);


main(int argc, char *argv[]){

  char signalList[CHARMAX];  //input name of signals file
  char ResultFileName[CHARMAX]; //File name to save rsults. If File exits append results, if no File is created
  int shotNumber;  // Input Shot number
  int Nsignals;    //Number the signal to use input+calculated

  signal *pSignal; //Pointer to struct signal, this struct contains all information of a signal

  FILE *fp;   //temporal pointer for file access
  int nlin;   //Counts the number of lines in file signallist. So the numebr of a raw signals
  char     txt[CHARMAX];   //Temporal array to hold path names array for input
  char   *ptr;   // temporal pointer to string
  char  title[53];
  char  units[11];

  float *SignalJET; //Para paso intermedio
  float *TimeJET; //Para paso intermedio


  int *t0;
  double *Maximums;
  double *Minimums;


  int i,j,t,z,k;  // for loops
  int error;  //Error indicator
  int ndata; //Number of datos of a signal

  double **pArray; //pointer to array
  double *pBuffer; //Pointer too buffer
  double *pBuffer1; //Pointer too buffer
  double *pBuffer2; //Pointer too buffer

  int  *ToNormalize; // pointer to array with actions about normalize proccess
                     //[12]={TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE};
  char **PathMx;  //Path to files M for model

  model *pModel;  //Pointer to struct model

  double tini;
  double torigin; //Initial time to resampling

  double dummy;
  double **Mdummy;

  int **rowptr;

  double R[5];   //Este hay que pasarlo tambien a dinamico
  double ModelParts[11][3];//Este hay que pasarlo tambien a dinamico
  double **ModelParts2;

  int NProcSignals;	//Number of processing signals =

  double ColumnModel[11];
  double D[3];	//tambien a dinÃ¡mico

  double Result_Model;
  int finish;


  //  double (*infft)[2];   //Buffer for FFT
  //  double (*outfft)[2];

  double infft[32][2];//ojo hay que pasarlo a dinmico
  double outfft[32][2];

  double resabs[32];

  int iWindow=0;
  int ndummy=0;
  int ndummy2=0;

	FILE *filelog;
	FILE *fileRes;


  /* Test of input parameters */
  if(argc != 4)
    {
      printf("\nProgram usage:\n\n");
      printf("> Main signalList shotNumber\n\n");
      printf("signalList: file with the JPF signal paths\n");
      printf("shotNumber: discharge to process\n");
      printf("Result File: File to append results \n");
      printf("\n\n");
      exit(0);
    }

  /* Input parameters to program variables */
  strcpy(signalList, argv[1]);
  shotNumber = atoi(argv[2]);
  strcpy(ResultFileName, argv[3]);


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
      //      printf("Valor de txt %s \n",txt);

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
    //    printf("Leido del struct %s \n",(pSignal+i)->name);
  }

  for (i=0;i<conf_Nsignals;i++){
    //    printf("Reading %s from shot %d\n",(pSignal+i)->name, shotNumber);
    fflush(stdout);

    getnwds_((pSignal+i)->name, &shotNumber, &ndata, &error, (long) strlen((pSignal+i)->name)); //Determine the space needed

    //    printf ("Debugg nData: %d error: %d\n", ndata, error);
    if ((!error) && (ndata > 0)){  //data are available
      if( ((pSignal+i)->pData = (double *) malloc(ndata*sizeof(double))) == (double *) 0 ){ //Allocate memory
              free(pSignal);
	      printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
	      exit(0);
       }
      if(((pSignal+i)->pTime = (double *) malloc(ndata*sizeof(double))) == (double *) 0){ //Allocate memory
	      free((pSignal+i)->pData);
              free(pSignal);
	      printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
	      exit(0);
       }

      if( (SignalJET = (float *) malloc(ndata*sizeof(float))) == (float *) 0 ){ //Allocate memory
              free(pSignal);
	      printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
	      exit(0);
       }
      if((TimeJET = (float *) malloc(ndata*sizeof(float))) == (float *) 0){ //Allocate memory
	      free((pSignal+i)->pData);
              free(pSignal);
	      printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
	      exit(0);
       }



       // Pre set last char of returned text to 0 for c string
       title [sizeof(title) - 1] = 0;
       units [sizeof(units) - 1] = 0;
       //  getdat_((pSignal+i)->name, &shotNumber, (pSignal+i)->pData, (pSignal+i)->pTime,&ndata , title, units, &error,(long) strlen((pSignal+i)->name), sizeof(title) - 1l, sizeof(units) - 1l);
       getdat_((pSignal+i)->name, &shotNumber, SignalJET, TimeJET ,&ndata , title, units, &error,(long) strlen((pSignal+i)->name), sizeof(title) - 1l, sizeof(units) - 1l);

       (pSignal+i)->nSamples= ndata;
       (pSignal+i)->Npoints= conf_Npoints;


       for (j=0;j<ndata; j++){
	 *((pSignal+i)->pData+j)= (double) *(SignalJET+j);
	 *((pSignal+i)->pTime+j)= (double) *(TimeJET+j);

       }


    	 txt[0]= NULL;
      	 sprintf(txt,"Original_%d.txt",i);
 /*     	 salvaOutput (((pSignal+i)->pData),txt,ndata);
      	 txt[0]= NULL;
      	 sprintf(txt,"O_Tiempos_%d.txt",i);
      	 salvaOutput (((pSignal+i)->pTime),txt,ndata);*/

      	  salvaResampling ( ((pSignal+i)->pData), ((pSignal+i)->pTime), txt, ndata);



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

  /*  for (i=1143;i<1155;i++){
    printf ("Signal t: %f valor %f \n",*((pSignal+0)->pTime+i),*((pSignal+0)->pData+i));

    }*/






  Maximums= (double *) malloc(sizeof(double)*Nsignals); //Alocate space for Maximums and Minimums
  Minimums= (double *) malloc(sizeof(double)*Nsignals); //Alocate space for Maximums and Minimums
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
   //      printf("Leido del struct %d \n",(pSignal+i)->Normalize);
   //  printf("VMax: %.10f  VMin:  %.10f \n",*(Maximums+i),*(Minimums+i));
  }

 free(Maximums); //free array values are copied in the struct
 free(Minimums);
 free(ToNormalize);



	filelog = fopen("W.txt", "w");
	fileRes = fopen( ResultFileName, "a");

  //Look for t0 for the signals, the reference is Ipla signal 0


  t0= (int *) malloc(sizeof(int)*Nsignals); //Alocate space for initial times in the signals

  *t0= IndexEvent((pSignal)->pData,(pSignal)->nSamples,conf_Threshold,0); //Look for the time when Ipla < Threshold

//  *t0= 998; //a Huevo para ajustar con Sebas
   printf("Indice de t0 %d \n",*t0);
   printf("Valor de t0 indice %d valor %f \n",*(t0),*(pSignal->pTime+(*t0)));


  // t0 for all signals, Only for raw signals

  for (i=1;i < conf_Nsignals;i++){
    *(t0+i)=  IndexEvent( (pSignal+i)->pTime,(pSignal+i)->nSamples,*((pSignal)->pTime+(*t0)),1)-1;
     printf("Valor de t0 indice %d valor %f \n",*(t0+i),*((pSignal+i)->pTime+(*(t0+i))) );
  }

  //make pointer array and allocate memory for memory windows in models
  printf("Size del double: %d \n",sizeof(double));
  for (i=0;i<Nsignals;i++){ //Allocate array's pointer and pointer to resamplig time
    (pSignal+i)->pM=  malloc(sizeof(double )*conf_NModels); //Allocate space for arrays pointers
    //(pSignal+i)->pM=  malloc(sizeof(double * )*conf_NModels); //Allocate space for arrays pointers  NO SERIA ESTA?
  //  printf("Puntero Para pSignal+%d ->pM: %p \n",i,(pSignal+i)->pM);
  }

  for (i=0;i<Nsignals;i++){ //Allocate memory buffers
    (pSignal+i)->pTimeR= malloc(sizeof(double)*conf_Npoints); //Buffer for resampling time
    for(j=0;j<conf_NModels;j++){
      *((pSignal+i)->pM+j)= (double *) malloc(conf_Npoints*sizeof(double));
      if(*((pSignal+i)->pM+j) == NULL ){ //Do so for easy code ready
	free(pSignal);
        printf("\n**** Not available memory for Model windows  %s, \n",(pSignal+i)->name );
        exit(0);
      }
    }
  }

  //Here the raw signal are read, and all buffers are allocated for all signal
  //Now the buffer are going to be filled with resampling data window.

  torigin= *((pSignal)->pTime+*(t0)); //Fix to Signal 1 origin

  tini= torigin; //Fix to Signal 1 origin

  for (j=conf_NModels-1;j>=0;j--){

	  for(z=0;z < conf_Npoints;z++){

		  for (i=0;i<conf_Nsignals;i++){ //Resampling for raw signals only
			     t=resampling2 ( *(t0+i), conf_Sampling, z, tini,*((pSignal+i)->pM+j) ,(pSignal+i));
		  }
		  tini += conf_Sampling;

	  }
	  for (i=0;i<conf_Nsignals;i++){ //Resampling for raw signals only
		  txt[0]= NULL;
		  sprintf(txt,"Resampleada_%s.txt",(pSignal+i)->name);

		  salvaResampling ( *((pSignal+i)->pM+j), (pSignal+i)->pTimeR, txt, (pSignal+i)->Npoints);

	  }

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
    txt[0]= NULL;
    sprintf(txt,"Signal_%d_%d.txt",7,j);
    salvaOutput (*((pSignal+7)->pM+j),txt,32);
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
    txt[0]= NULL;
    sprintf(txt,"Signal_%d_%d.txt",7,j);
    salvaOutput (*((pSignal+7)->pM+j),txt,32);
  }



  //Signal 8. Done

  //Signal 9. difference between signals 6 and 7. The signal is number 8 (9-1)

  for (j=conf_NModels-1;j>0;j--){  //two parts loop for use the last loop M3 to generate the timeR
     // printf ("pM %d : ",(pSignal+i)->pM+j);
    pBuffer= *((pSignal+8)->pM+j);  //Used for easy code reading
    pBuffer1= *((pSignal+5)->pM+j);
    pBuffer2= *((pSignal+6)->pM+j);
    for (i=0;i<conf_Npoints;i++){
      *(pBuffer+i)= normalize ((pSignal+8)->Max, (pSignal+8)->Min, (*(pBuffer2+i) < 1 ? 1 :*(pBuffer2+i)) / (*(pBuffer1+i) < 1000 ? 1000 :*(pBuffer1+i))); //Ojo 1000  y 1 a configuracion
    }
    txt[0]= NULL;
    sprintf(txt,"Signal_%d_%d.txt",8,j);
    salvaOutput (*((pSignal+8)->pM+j),txt,32);
  }


  for (j=0;j>=0;j--){                           //genarate M3 timeR and data
     // printf ("pM %d : ",(pSignal+i)->pM+j);
    pBuffer= *((pSignal+8)->pM+j);  //Used for easy code reading
    pBuffer1= *((pSignal+5)->pM+j);
    pBuffer2= *((pSignal+6)->pM+j);

    for (i=0;i<conf_Npoints;i++){
      *(pBuffer+i)= normalize ((pSignal+8)->Max, (pSignal+8)->Min, (*(pBuffer2+i) < 1 ? 1 :*(pBuffer2+i)) / (*(pBuffer1+i) < 1000 ? 1000 :*(pBuffer1+i))); //Ojo 1000  y 1 a configuracion
      *((pSignal+8)->pTimeR+i)= *((pSignal+5)->pTimeR+i);
    }
    txt[0]= NULL;
    sprintf(txt,"Signal_%d_%d.txt",8,j);
    salvaOutput (*((pSignal+8)->pM+j),txt,32);
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

     (pModel+i)->alfa = malloc(sizeof(double )*(pModel+i)->nvectors); //Allocate space for arrays
     if (Mdummy == NULL){
        free(pSignal);
        printf ("Not space for alfa data array \n");
        exit(0);
     }

     Mdummy = malloc(sizeof(double *)*(pModel+i)->nvectors); //Allocate space for arrays of data

     if (Mdummy == NULL){
       free(pSignal);
       printf ("Not space for Mmodel data array \n");
       exit(0);
     }

     for (j=0;j<(pModel+i)->nvectors;j++){
       Mdummy[j]= malloc(sizeof(double)*(pModel+i)->coef_vector); //Allocate space for arrays
       if (Mdummy[j] == NULL){
	 free(pSignal);
	 printf ("Not space for Mmodel data array \n");
	 exit(0);
       }

     }

     (pModel+i)->data= Mdummy;

     M_values (*(PathMx+i), pModel+i); //Load model firts M1 M2 M3


     //    printf("\n Vectot m1 \n");
     //for(t=0;t<11;t++)
     //printf("%.10f \t",*(Mdummy+t));
   }





   Read_M (conf_PathR, R);

   printf ("The model has been read \n");


   //here the model is read

   //calculamos los valores del medelo de las tres ventanas
   //esto hay que meterlo dentro de la estructira de punteros

   ModelParts2= malloc(sizeof(double * )* (Nsignals*2));
   if (ModelParts2 == NULL){
      free(pSignal);
      printf ("Not space for ModelParts data array \n");
      exit(0);
   }
   for (i=0;i<(Nsignals*2);i++){
	   *(ModelParts2+i)= malloc(sizeof(double ) * conf_NModels);
	   if (*(ModelParts2+i) == NULL){
	      free(pSignal);
	      printf ("Not space for ModelParts data array \n");
	      exit(0);
	   }
   }

   //(pSignal+i)->pM=  malloc(sizeof(double * )*conf_NModels); //Allocate space for arrays pointers  NO SERIA ESTA?



//Process Ipla


/*   for (i=0;i<conf_NModels;i++){ //mean ipla
     ModelParts[0][i]= Mean (*((pSignal)->pM+i), conf_Npoints);
     *(*(ModelParts2)+i)= Mean (*((pSignal)->pM+i), conf_Npoints);
   }

     salvaOutput (pBuffer,"FFTinIPLA.txt",32);

*/



  //Process loca

   for(j=0;j<Nsignals*2;j+2){

	   for (i=0;i<conf_NModels;i++){ //mean values
	        //ModelParts[2][i]= Mean (*((pSignal+1)->pM+i), conf_Npoints);
	        *(*(ModelParts2+j)+i)= Mean (*((pSignal)->pM+i), conf_Npoints);
	   }


	   for (i=0;i<conf_NModels;i++){ //des loca
	        pBuffer= *((pSignal+1)->pM+i);
	        for (t=0;t<conf_Npoints;t++){
	          infft[t][0]=*(pBuffer+t); infft[t][1]=0;
	        }

	        fft(conf_Npoints, infft, outfft);

	        Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

	        //ModelParts[3][i]= Desv (resabs, conf_Npoints/2);
	        *(*(ModelParts2+(j+1))+i)= Desv (resabs, conf_Npoints/2);

	      }



   }


  //para probrar distancia
  /* ColumnModel[0]=  0.666;
  ColumnModel[1]=  0.00524;
  ColumnModel[2]=  0.01588;
  ColumnModel[3]=  0.00301;
  ColumnModel[4]=  0.001;
  ColumnModel[5]=  0.41767;
  ColumnModel[6]=  0.00218;
  ColumnModel[7]=  0.50595;
  ColumnModel[8]=  0.00033;
  ColumnModel[9]=  0;
  ColumnModel[10]= 0.00002;*/


  D[0]=  distance (ColumnModel, pModel);
  // printf ("Distancia: %.10f \n",D[0]);

  // dummy= prod_vect (ColumnModel, *pModel->data, pModel->bias, 11);

  //printf ("Distancia dummy : %.10f \n",dummy);



  for (i=0;i<11;i++){
    ColumnModel[i]= ModelParts[i][1];
  }
  D[1]=  distance (ColumnModel, pModel+1);

  for (i=0;i<11;i++){
    ColumnModel[i]= ModelParts[i][2];
  }

  D[2]=  distance (ColumnModel, pModel+2);
  //  printf ("Distancia: %.10f \n",D[0]);
  //  printf ("Distancia: %.10f \n",D[1]);
  // printf ("Distancia: %.10f \n",D[2]);

  //  printf("tR: %f \t %f  \t  %f  \t %f \n",*((pSignal)->pTimeR+31),*((pSignal+2)->pTimeR+31),*((pSignal+3)->pTimeR+31),*((pSignal+7)->pTimeR+31));

  fprintf(filelog,"tiempo \t S1 \t\t S2 \t S3 \t S4 \t\t S5 \t\t S6 \t S7 \t\t S8 \tS9 \t\t S10 \t S11 \t D3 \t   D2 \t    D1 \t      R\n");


  // printf ("\n Reading distance \n");
  //  fprintf(filelog,"vector D: %.10f \t %.10f \t %.10f \n",D[2],D[1],D[0]);

  Result_Model= (double) ((D[2]*R[2] + D[1]*R[3] + D[0]*R[4] +R[1])/R[0]);
  //printf("Resultado D2: %.10f  D1: %.10f  D0: %.10f \n",D[2],D[1],D[0]);
  // fprintf(filelog,"%.10f \n",Result_Model);

  //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][2],ModelParts[1][2],ModelParts[2][2],ModelParts[3][2],ModelParts[4][2],ModelParts[5][2],ModelParts[6][2],ModelParts[7][2],ModelParts[8][2],ModelParts[9][2],ModelParts[10][2],D[2],D[1],D[0],Result_Model);
  //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][1],ModelParts[1][1],ModelParts[2][1],ModelParts[3][1],ModelParts[4][1],ModelParts[5][1],ModelParts[6][1],ModelParts[7][1],ModelParts[8][1],ModelParts[9][1],ModelParts[10][1],D[2],D[1],D[0],Result_Model);
  fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][0],ModelParts[1][0],ModelParts[2][0],ModelParts[3][0],ModelParts[4][0],ModelParts[5][0],ModelParts[6][0],ModelParts[7][0],ModelParts[8][0],ModelParts[9][0],ModelParts[10][0],D[2],D[1],D[0],Result_Model);




//Evaluar result si es activo escribir alarma y exit
// si no bucle





  //    N_vectors("p.txt" , &dmmy, &dmmy1);

  //   printf ("resultados %d %d \n",dmmy,dmmy1);

   /**************************************
   // LOOP
   ******************************/




   finish= FALSE;
k=3;
   do{
	 ndummy=iwindow1(&iWindow,conf_NModels);
     j=bufferFree(&iWindow,conf_NModels);

	  for(z=0;z < conf_Npoints;z++){

		  for (i=0;i<conf_Nsignals;i++){ //Resampling for raw signals only
			  //if (i==0)    printf("tiempos indice t: %d \n", *(t0+i));
			  //tini= *((pSignal+i)->pTimeR+31); //the time always is pTimeR
			  //t=resampling ( *(t0+i), conf_Sampling,*((pSignal+ndummy)->pTimeR+31) ,*((pSignal+i)->pM+j) ,(pSignal+i));
			  //t=resampling ( *(t0+i), conf_Sampling,tini ,*((pSignal+i)->pM+j) ,(pSignal+i));
			 // *(t0+i)+=t;

			  /*       txt[0]= NULL;
			   sprintf(txt,"Signal_%d_%d.txt",i,k);
			   salvaOutput (*((pSignal+i)->pM+j),txt,32);
			   txt[0]= NULL;
			   sprintf(txt,"Tiempos_%d_%d.txt",i,k);
			   salvaOutput (((pSignal+i)->pTimeR),txt,32);
			   */
			   t=resampling2 ( *(t0+i), conf_Sampling, z, tini,*((pSignal+i)->pM+j) ,(pSignal+i));


			  /*
			  //       if (i==0) printf("tiempos indice t retorno : %d \n",t);
			  if (t==0){
				  printf("\n ***********  Signal end reach No disruption **************\n");
				  fprintf (fileRes,"%d \t %s \t %.3f \n",shotNumber,"-1",*((pSignal)->pTimeR+31));
				  finish= TRUE;
				  break;
			  }
			  */
		  }
		  tini += conf_Sampling;

	  }
	  printf("%.3f \n",tini);
	  for (i=0;i<conf_Nsignals;i++){ //Resampling for raw signals only
		  txt[0]= NULL;
		  sprintf(txt,"Resampleada_%s.txt",(pSignal+i)->name);

		  salvaResampling ( *((pSignal+i)->pM+j), (pSignal+i)->pTimeR, txt, (pSignal+i)->Npoints);

	  }

	  k++;
	  if (k>53)
		  k=105;

     //Here the raw signal are reading and ready to be used

     //Now start the caculate signals. Signal number 8 (from number 3) the difference signal
     //Signal 9 (from 6 and 7). Signal 6 is adjust to a minimum value of 1000 and then is divided
     //for signal 7.

     //Siganl 8. Index=7 (8 - 1)
     //Make the new signal x[n]=x[n+1]-x[n], the last point is the previous point

     pBuffer= *((pSignal+7)->pM+j);  //Used for easy code reading
     pBuffer1= *((pSignal+2)->pM+j);

     for (i=0;i<conf_Npoints-1;i++){
       *(pBuffer+i)= normalize ((pSignal+7)->Max, (pSignal+7)->Min, *(pBuffer1+i+1)-*(pBuffer1+i));
       *((pSignal+7)->pTimeR+i)= *((pSignal+2)->pTimeR+i);
     }
     *(pBuffer+i)= *(pBuffer+i-1);
     *((pSignal+7)->pTimeR+i)= *((pSignal+7)->pTimeR+i-1);



     //Signal 8. Done

     //Signal 9. difference between signals 6 and 7. The signal is number 8 (9-1)

     // printf ("pM %d : ",(pSignal+i)->pM+j);
     pBuffer= *((pSignal+8)->pM+j);  //Used for easy code reading
     pBuffer1= *((pSignal+5)->pM+j);
     pBuffer2= *((pSignal+6)->pM+j);

     for (i=0;i<conf_Npoints;i++){
       *(pBuffer+i)= normalize ((pSignal+8)->Max, (pSignal+8)->Min, (*(pBuffer2+i) < 1 ? 1 :*(pBuffer2+i)) / (*(pBuffer1+i) < 1000 ? 1000 :*(pBuffer1+i))); //Ojo 1000  y 1 a configuracion
       *((pSignal+8)->pTimeR+i)= *((pSignal+5)->pTimeR+i);
     }


/*		   txt[0]= NULL;
		   sprintf(txt,"Signal_8_%d.txt",k,32);
		   salvaOutput (*((pSignal+8)->pM+j),txt,32);
		   txt[0]= NULL;
		   sprintf(txt,"Tiempos_8_%d.txt",k,32);
		   salvaOutput (((pSignal+8)->pTimeR),txt,32);
		   k++;
*/



   //calculamos los valores del medelo de las tres ventanas
   //esto hay que meterlo dentro de la estructira de punteros


   //Process Ipla

   ModelParts[0][j]= Mean (*((pSignal)->pM+j), conf_Npoints);

   //des ipla
   pBuffer= *((pSignal)->pM+j);
   for (t=0;t<conf_Npoints;t++){
     infft[t][0]=*(pBuffer+t); infft[t][1]=0;
   }

   fft(conf_Npoints, infft, outfft);

   Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

   ModelParts[1][j]= Desv (resabs, conf_Npoints/2);

//Process loca


   ModelParts[2][j]= Mean (*((pSignal+1)->pM+j), conf_Npoints);

 //des loca
   pBuffer= *((pSignal+1)->pM+j);
   for (t=0;t<conf_Npoints;t++){
     infft[t][0]=*(pBuffer+t); infft[t][1]=0;
   }

   fft(conf_Npoints, infft, outfft);

   Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

   ModelParts[3][j]= Desv (resabs, conf_Npoints/2);



//Process density

 //mean density
   ModelParts[4][j]= Mean (*((pSignal+3)->pM+j), conf_Npoints);

//Process Desv. Energy

   ModelParts[5][j]= Mean (*((pSignal+4)->pM+j), conf_Npoints);

 //des Desv. Energy
   pBuffer= *((pSignal+4)->pM+j);
   for (t=0;t<conf_Npoints;t++){
     infft[t][0]=*(pBuffer+t); infft[t][1]=0;
   }

   fft(conf_Npoints, infft, outfft);

   Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

   ModelParts[6][j]= Desv (resabs, conf_Npoints/2);


//Process Der. Inductancia

   ModelParts[7][j]= Mean (*((pSignal+7)->pM+j), conf_Npoints);

//des der inductacia
   pBuffer= *((pSignal+7)->pM+j);
   for (t=0;t<conf_Npoints;t++){
     infft[t][0]=*(pBuffer+t); infft[t][1]=0;
   }

   fft(conf_Npoints, infft, outfft);

   Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

   ModelParts[8][j]= Desv (resabs, conf_Npoints/2);




//Process Cociente

 //mean quotien
   ModelParts[9][j]= Mean (*((pSignal+8)->pM+j), conf_Npoints);


 //des der inductacia
   pBuffer= *((pSignal+8)->pM+j);
   for (t=0;t<conf_Npoints;t++){
     infft[t][0]=*(pBuffer+t); infft[t][1]=0;
   }

     fft(conf_Npoints, infft, outfft);

     Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

     ModelParts[10][j]= Desv (resabs, conf_Npoints/2);



     //  printf("\n Valores signal: %d valor:  %.10f ", 0, ModelParts[0][iWindow]);
     //  fprintf(filelog,"%.10f \n",ModelParts[8][iWindow]);


  for (i=0;i<11;i++){
    ColumnModel[i]= ModelParts[i][iWindow];
  }
  //  D[0]=  distance (ColumnModel, pModel+iWindow);
  D[0]=  distance (ColumnModel, pModel); //El modelo siempre es el mismo que D es decir 0 en este caso

  ndummy=iwindow_1(&iWindow,conf_NModels);
  for (i=0;i<11;i++){
    ColumnModel[i]= ModelParts[i][ndummy];
  }
  // D[2]=  distance (ColumnModel, pModel+ndummy);
  D[2]=  distance (ColumnModel, pModel+2); //El modelo siempre es el mismo que D es decir 2 en este caso

  ndummy2=ndummy;

  ndummy= iwindow_1(&ndummy,conf_NModels);
  for (i=0;i<11;i++){
    ColumnModel[i]= ModelParts[i][ndummy];
  }
  // D[1]=  distance (ColumnModel, pModel+ndummy);
  D[1]=  distance (ColumnModel, pModel+1); //El modelo siempre es el mismo que D es decir 1 en este caso





  //  fprintf(filelog,"vector D: %.10f \t %.10f \t %.10f \n",D[2],D[1],D[0]);

   Result_Model= (double) ((D[2]*R[2] + D[1]*R[3] + D[0]*R[4] +R[1])/R[0]);
   //   Result_Model= D[2]*R[3] + D[1]*R[2] + D[0]*R[1] +R[0];
   //   printf("Result_Model: %.5f \n",Result_Model);
   // fprintf(filelog,"%.10f \n",Result_Model);

   //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][ndummy2],ModelParts[1][ndummy2],ModelParts[2][ndummy2],ModelParts[3][ndummy2],ModelParts[4][ndummy2],ModelParts[5][ndummy2],ModelParts[6][ndummy2],ModelParts[7][ndummy2],ModelParts[8][ndummy2],ModelParts[9][ndummy2],ModelParts[10][ndummy2],D[2],D[1],D[0],Result_Model);
   //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][ndummy],ModelParts[1][ndummy],ModelParts[2][ndummy],ModelParts[3][ndummy],ModelParts[4][ndummy],ModelParts[5][ndummy],ModelParts[6][ndummy],ModelParts[7][ndummy],ModelParts[8][ndummy],ModelParts[9][ndummy],ModelParts[10][ndummy],D[2],D[1],D[0],Result_Model);
   fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][iWindow],ModelParts[1][iWindow],ModelParts[2][iWindow],ModelParts[3][iWindow],ModelParts[4][iWindow],ModelParts[5][iWindow],ModelParts[6][iWindow],ModelParts[7][iWindow],ModelParts[8][iWindow],ModelParts[9][iWindow],ModelParts[10][iWindow],D[2],D[1],D[0],Result_Model);

   printf ("ventana iwindows: %d  ndummy: %d  ndummy2:  %d  \n", iWindow, ndummy, ndummy2);



   if (Result_Model > 0){
     printf("\n ********* Disruption at  t: %f       *********** \n",*((pSignal)->pTimeR+31));
     fprintf (fileRes,"%d \t %s \t %.3f \n",shotNumber,"+1",*((pSignal)->pTimeR+31));
     finish=TRUE;
   }




   }while(!finish );


	fclose(filelog);
	fclose(fileRes);







}





void PRINTSIGNAL( signal *entrada){
  int i;

  printf ("Valores comunes Max: %f  Min: %f Pts: %d \n",entrada->Max, entrada->Min, entrada->Npoints);
  printf ("ind \t Tiempo \t\t valor M3 \t valor M2 \t valor M1 \n");
  for (i=0;i<entrada->Npoints;i++){
    printf ("%d \t tr: %f  \t  %f  \t  %f \t  %f \n", i, *(entrada->pTimeR+i),*(entrada->pM+2+i),*(entrada->pM+1+i),*(entrada->pM+i));
  }

}
