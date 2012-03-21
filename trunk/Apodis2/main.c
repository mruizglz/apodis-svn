
 /*

Program implemented by Juan Manuel Lopez, Universidad Politecnica de Madrid
july-Augoust 2011
Implementation of Apodis 2

JET signal reading method provide by Jesus Vega (CIEMAT)

Se ha generado un branch para estabilizar antes de hacer los modelos de 8 y 12

*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iso646.h>
#include <ppf.h>
#include "tipos.h"
#include "getfix/getfix.h"
#include "./tools/tools.h"
#include "./configuracion/config.h"
#include "./fft/fft.h"



void PRINTSIGNAL( signal *entrada);

//#define DEBUGLEVEL1	1
//#define DEBUGLEVEL2	2
//#define DEBUGSINCRO	3
//#define DEBUGWINDOW	4
//#define DEBUGWNAME	5


//prueba solo verificacion funciona el svn

main (int argc, char *argv[]){

  char signalList[CHARMAX];  //input name of signals file
  char ResultFileName[CHARMAX]; //File name to save rsults. If File exits append results, if no File is created
  int shotNumber;  // Input Shot number
  int Nsignals;    //Number the signal to use input+calculated

  char procesada[CHARMAX];

  signal *pSignal; //Pointer to struct signal, this struct contains all information of a signal

  FILE *fp;   //temporal pointer for file access
  int nlin;   //Counts the number of lines in file signallist. So the numebr of a raw signals
  char     txt[CHARMAX];   //Temporal array to hold path names array for input
  char   *ptr;   // temporal pointer to string
  char  title[53];
  char  units[11];

  float *SignalJET; //Para paso intermedio
  float *TimeJET; //Para paso intermedio


  int *t0;			//Pointer to init times, cross threshold array
  double *Maximums;	//Pointer to maximums array
  double *Minimums;	//Pointer to minimums array


  int i,j,t,z,k;  // for loops
  int error;  //Error indicator
  int ndata; //Number of datos of a signal

//  double **pArray; //pointer to array
  double *pBuffer; //Pointer to buffer
//  double *pBuffer1; //Pointer to buffer
//  double *pBuffer2; //Pointer to buffer

  int  *ToNormalize; // pointer to array with actions about normalize proccess
                     //[12]={TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE};
  char **PathMx;  //Path to files M for model

  model *pModel;  //Pointer to struct model

  double tini;
  double torigin; //Initial time to resampling

  double dummy;
  double **Mdummy;

//  int **rowptr;

  int Ncoefficients;

  double R[5];   //Este hay que pasarlo tambien a dinamico
//  double ModelParts[11][3];//Este hay que pasarlo tambien a dinamico
  double **ModelParts2;

  //int NProcSignals;	//Number of processing signals =

  //double ColumnModel[11];
  //double D[3];	//tambien a dinÃ¡mico
  double *ColumnModel;
  double *D;


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

  double finaltime;
  int tcurrent;

  //Variables for PPF processing

  char *pPwd = "";
  char *pUser = "ppfread";
  int status;
  int lun = 0;
  int seqNumber = 0;
  int shotppf = 0;
  int nx,nt,systat,ustat;
  char ppfsignal[128], dda[128];
  char *pslash;
  char frmt[8],comm[32],ihdat[60];
  int irdat[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int iwdat[13];
  float *pX0;



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

  #ifdef DEBUGLEVEL1
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
  #endif

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

   if ((pSignal+i)->name[0] == '*' ){ //PPF signal's name start by *

    	strcpy((pSignal+i)->name,(pSignal+i)->name+1); //delete *

    	PPFPWD(pUser, pPwd, (unsigned int) strlen(pUser), (unsigned int) strlen(pPwd));
    	PPFPOK(&status);
    	if(status){
    	      printf("Error in ppfpwd(). Status: %d\n\nEnd of program\n", status);
    	      exit(0);
    	}
    	PPFGO(&shotppf, &seqNumber, &lun, &status);
    	if(status){
    	      printf("Error in ppfgo(). Status: %d\n\nEnd of program\n", status);
    	      exit(0);
    	}

    	strcpy(dda,(pSignal+i)->name);
    	pslash= strstr(dda,"/");
    	strcpy(ppfsignal,pslash+1);
    	*pslash=NULL;

    	PPFDTI(&shotNumber, &seqNumber, dda, ppfsignal, &nx, &nt, frmt, units, comm, &systat, &ustat, &status,
    		 (unsigned int) sizeof(dda),  (unsigned int) sizeof(signal),  (unsigned int) sizeof(frmt),
    		 (unsigned int) sizeof(units),  (unsigned int) sizeof(comm));
        if(status){
    	      printf("Error in ppdti(). Status: %d\n\nEnd of program\n", status);
    	      exit(0);
    	}
        irdat[1] = nx;
        irdat[3] = nt;
        irdat[5] = irdat[1] - irdat[0];
        irdat[6] = irdat[3] - irdat[2];
        irdat[4] = irdat[5]*irdat[6];

		if ((!status) && (irdat[4] > 0)){  //data are available
			if( ((pSignal+i)->pData = (double *) malloc(irdat[4]*sizeof(double))) == (double *) 0 ){ //Allocate memory
				free(pSignal);
				printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
				exit(0);
			}
			if(((pSignal+i)->pTime = (double *) malloc(irdat[6]*sizeof(double))) == (double *) 0){ //Allocate memory
				free((pSignal+i)->pData);
				free(pSignal);
				printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
				exit(0);
			}

			if( (SignalJET = (float *) malloc(irdat[4]*sizeof(float))) == (float *) 0 ){ //Allocate memory
				free(pSignal);
				printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
				exit(0);
			}
			if((TimeJET = (float *) malloc(irdat[6]*sizeof(float))) == (float *) 0){ //Allocate memory
				free((pSignal+i)->pData);
				free(pSignal);
				printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+i)->name, shotNumber);
				exit(0);
			}

			if((pX0 = (float *) malloc(irdat[5]*sizeof(float))) == (float *) 0)	{
				  free((pSignal+i)->pData);
				  free(pSignal);
			      printf("Error in malloc(X)\n\nEnd of program\n");
			      exit(0);
			}

			/*   memset(ihdat, 0, sizeof(ihdat)); */
			  PPFGET(&shotNumber, dda, ppfsignal, irdat, ihdat, iwdat, SignalJET, pX0, TimeJET, &status, (unsigned int) sizeof(dda),
				 (unsigned int) sizeof(signal), (unsigned int) sizeof(ihdat));
			  if(status)
			    {
				  free((pSignal+i)->pData);
				  free(pSignal);
			      printf("Error in ppfget(). Status: %d\n\nEnd of program\n", status);
			      exit(0);
			    }


			(pSignal+i)->nSamples= irdat[4];
			(pSignal+i)->Npoints= conf_Npoints;


			for (j=0;j<irdat[4]; j++){
				*((pSignal+i)->pData+j)= (double) *(SignalJET+j);
				*((pSignal+i)->pTime+j)= (double) *(TimeJET+j);

			}

			#ifdef DEBUGLEVEL1
				txt[0]= NULL;
				sprintf(txt,"Original_%d.txt",i);


				/*     	 salvaOutput (((pSignal+i)->pData),txt,ndata);
				 txt[0]= NULL;
				 sprintf(txt,"O_Tiempos_%d.txt",i);
				 salvaOutput (((pSignal+i)->pTime),txt,ndata);*/

				salvaResampling ( ((pSignal+i)->pData), ((pSignal+i)->pTime), txt, irdat[4]);
			#endif


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
    else {
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

			#ifdef DEBUGLEVEL1
				txt[0]= NULL;
				sprintf(txt,"Original_%d.txt",i);


				/*     	 salvaOutput (((pSignal+i)->pData),txt,ndata);
				 txt[0]= NULL;
				 sprintf(txt,"O_Tiempos_%d.txt",i);
				 salvaOutput (((pSignal+i)->pTime),txt,ndata);*/

				salvaResampling ( ((pSignal+i)->pData), ((pSignal+i)->pTime), txt, ndata);
			#endif


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


  } //end for

  //Model 12. Signal 7 is Poloidal beta; Signal 8 Vertical position. Signal 9 derivate of 2 (inductance) Signal 10 derivate 7 (poloidal beta) Signal 11 derivate 8 (Vertical Position)


	strcpy((pSignal+9)->name,"Der. Inductance"); //Signal 9 Der. Inductance
	if( ((pSignal+9)->pData = (double *) malloc((pSignal+2)->nSamples*sizeof(double))) == (double *) 0 ){ //Allocate memory
		free(pSignal);
		printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+9)->name, shotNumber);
		exit(0);
	}
	if(((pSignal+9)->pTime = (double *) malloc((pSignal+2)->nSamples*sizeof(double))) == (double *) 0){ //Allocate memory
		free((pSignal+9)->pData);
		free(pSignal);
		printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+9)->name, shotNumber);
		exit(0);
	}
	(pSignal+9)->nSamples= (pSignal+2)->nSamples;
	(pSignal+9)->Npoints= conf_Npoints;

 for (j=0;j<(pSignal+2)->nSamples-1;j++){
//		 dummy= (double) (*((pSignal+4)->pData+j+1) - *((pSignal+4)->pData+j) );
		*((pSignal+9)->pData+j)= (double) (*((pSignal+2)->pData+j+1) - *((pSignal+2)->pData+j) );
		*((pSignal+9)->pTime+j)= (double) (*((pSignal+2)->pTime+j+1) - *((pSignal+2)->pTime+j) );
  }
  printf("Signal 10 is derivate  \n");


	strcpy((pSignal+10)->name,"Der. Beta"); //Signal 10 Der. Poloidal beta
	if( ((pSignal+10)->pData = (double *) malloc((pSignal+7)->nSamples*sizeof(double))) == (double *) 0 ){ //Allocate memory
		free(pSignal);
		printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+10)->name, shotNumber);
		exit(0);
	}
	if(((pSignal+10)->pTime = (double *) malloc((pSignal+7)->nSamples*sizeof(double))) == (double *) 0){ //Allocate memory
		free((pSignal+10)->pData);
		free(pSignal);
		printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+10)->name, shotNumber);
		exit(0);
	}
	(pSignal+10)->nSamples= (pSignal+7)->nSamples;
	(pSignal+10)->Npoints= conf_Npoints;

for (j=0;j<(pSignal+7)->nSamples-1;j++){
//		 dummy= (double) (*((pSignal+4)->pData+j+1) - *((pSignal+4)->pData+j) );
		*((pSignal+10)->pData+j)= (double) (*((pSignal+7)->pData+j+1) - *((pSignal+7)->pData+j) );
		*((pSignal+10)->pTime+j)= (double) (*((pSignal+7)->pTime+j+1) - *((pSignal+7)->pTime+j) );
}
printf("Signal 11 is derivate  \n");

strcpy((pSignal+11)->name,"Der. Vertical Pos."); //Signal 11 Der. Vertical Pos
if( ((pSignal+11)->pData = (double *) malloc((pSignal+8)->nSamples*sizeof(double))) == (double *) 0 ){ //Allocate memory
	free(pSignal);
	printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+11)->name, shotNumber);
	exit(0);
}
if(((pSignal+11)->pTime = (double *) malloc((pSignal+8)->nSamples*sizeof(double))) == (double *) 0){ //Allocate memory
	free((pSignal+11)->pData);
	free(pSignal);
	printf("\n**** Not available memory for signal %s, shot %d\n\n",(pSignal+11)->name, shotNumber);
	exit(0);
}
(pSignal+11)->nSamples= (pSignal+8)->nSamples;
(pSignal+11)->Npoints= conf_Npoints;

for (j=0;j<(pSignal+8)->nSamples-1;j++){
//		 dummy= (double) (*((pSignal+4)->pData+j+1) - *((pSignal+4)->pData+j) );
	*((pSignal+11)->pData+j)= (double) (*((pSignal+8)->pData+j+1) - *((pSignal+8)->pData+j) );
	*((pSignal+11)->pTime+j)= (double) (*((pSignal+8)->pTime+j+1) - *((pSignal+8)->pTime+j) );
}
printf("Signal 12 is derivate  \n");





  double *pData;        //Pointer to raw data
  double *pTime;        //pointer to time for raw data signal
  double *pTimeR;       //pointer to resampling times for nearest time window
  double **pM;         //pointer to array of pointer of model windows
  int Normalize;        //the signal has to be Normalize

  //At this point the signal have been read form the JET BBDD

  //Handle exceptions

  for (i=0;i < (pSignal+5)->nSamples;i++){
	  *((pSignal+5)->pData+i)= *((pSignal+5)->pData+i) < 1000 ? 1000 : *((pSignal+5)->pData+i);
  }
  for (i=0;i<(pSignal+6)->nSamples;i++){
	  *((pSignal+6)->pData+i)= *((pSignal+6)->pData+i) < 1 ? 1 : *((pSignal+6)->pData+i);
  }

  finaltime= *(pSignal->pTime+((pSignal->nSamples)-5)); //final time for Ipla minus 5 samples for security

  /*  for (i=1143;i<1155;i++){
    printf ("Signal t: %f valor %f \n",*((pSignal+0)->pTime+i),*((pSignal+0)->pData+i));

    }*/






  Maximums= (double *) malloc(sizeof(double)*Nsignals); //Alocate space for Maximums and Minimums
  Minimums= (double *) malloc(sizeof(double)*Nsignals); //Alocate space for Maximums and Minimums
  ToNormalize= (int *) malloc(sizeof(int)*Nsignals); //Alocate space for normaliza Flag

  // Reading Maximums and Minimums values

  // Check Maximums, Minimum and Normalize sizes

  if((ReadFloatTxt(conf_PathMax, Maximums) != Nsignals) or (ReadFloatTxt(conf_PathMin, Minimums) !=  Nsignals or (ReadNormalizeTxt(conf_PathNormalize, ToNormalize) != Nsignals))){
    printf ("Maximum, Minimus, or Normalize  Files does not match with signal numbers \n");
    exit (0);
  }

 for(i=0;i<Nsignals;i++){	//Copy Maximums, minimums and normaliza to signaal struct
   (pSignal+i)->Max= *(Maximums+i);
   (pSignal+i)->Min= *(Minimums+i);
   (pSignal+i)->Normalize= *(ToNormalize+i);
   //      printf("Leido del struct %d \n",(pSignal+i)->Normalize);
   //  printf("VMax: %.10f  VMin:  %.10f \n",*(Maximums+i),*(Minimums+i));
  }

 free(Maximums); //free arrays, values are copied in the struct
 free(Minimums);
 free(ToNormalize);

	sprintf(procesada,"W.txt"); //Reuse procesada buffer
	#ifdef DEBUGWNAME
		sprintf(procesada,"W%d.txt",shotNumber); //Reuse procesada buffer
	#endif

	filelog = fopen(procesada, "w");
	fileRes = fopen( ResultFileName, "a");

  //Look for t0 for the signals, the reference is Ipla signal 0


  t0= (int *) malloc(sizeof(int)*Nsignals); //Alocate space for initial times in the signals

  *t0= IndexEvent((pSignal)->pData,(pSignal)->nSamples,conf_Threshold,0); //Look for the time when Ipla < Threshold

  //*t0= 1144; //a Huevo para ajustar con Sebas
   //printf("Indice de t0 %d \n",*t0);
   printf("Valor de t0 indice %d valor %f \n",*(t0),*(pSignal->pTime+(*t0)));


  // t0 for all signals, Only for raw signals

  for (i=1;i < conf_Nsignals;i++){
    *(t0+i)=  IndexEvent( (pSignal+i)->pTime,(pSignal+i)->nSamples,*((pSignal)->pTime+(*t0)),1)-1;
     //printf("Valor de t0 indice %d valor %f \n",*(t0+i),*((pSignal+i)->pTime+(*(t0+i))) );
  }

  //make pointer array and allocate memory for memory windows in models
  //printf("Size del double: %d \n",sizeof(double));
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
  //tini= 42.330061; //Fix to Signal 1 origin

  tcurrent= IndexEvent(((pSignal)->pData+(*t0)+50),((pSignal)->nSamples)-((*t0)+50),conf_Threshold,1); //Look for the time when Ipla < Threshold
  finaltime= *(pSignal->pTime+(tcurrent+(*t0)+50));

	#ifdef DEBUGSINCRO
	  //We are to looking for the initial time in sebas's files

	  sprintf(procesada,"./procesadas/DES_%d_01_proc.txt",shotNumber);

	  tini= tsincro (procesada, tini);

	  printf("Sincronize signal at : %f \n",tini);
	#endif


  for (j=conf_NModels-1;j>=0;j--){	//Fill

	  for(z=0;z < conf_Npoints;z++){

		  for (i=0;i<Nsignals;i++){ //Resampling for raw signals only  =====>  se pasa de conf_Nsignals a Nsignals
			     t=resampling2 ( *(t0+i), conf_Sampling, z, tini,*((pSignal+i)->pM+j) ,(pSignal+i));
		  }
		  tini += conf_Sampling;

	  }
	  for (i=0;i<Nsignals;i++){ //Resampling for raw signals only  =====>  se pasa de conf_Nsignals a Nsignals


		  #ifdef DEBUGLEVEL1
			  txt[0]= NULL;
			  sprintf(txt,"Resampleada_%s.txt",(pSignal+i)->name);

			  salvaResampling ( *((pSignal+i)->pM+j), (pSignal+i)->pTimeR, txt, (pSignal+i)->Npoints);
          #endif
	  }

  }





  //Here the raw signal are reading and ready to be used




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

	 #ifdef DEBUGLEVEL1
		 printf("Tamanos Modelo: %d vectores: %d coeficientes: %d \n",i,(pModel+i)->nvectors , (pModel+i)->coef_vector);
	 #endif

     (pModel+i)->alfa = malloc(sizeof(double )*(pModel+i)->nvectors); //Allocate space for arrays
     if ((pModel+i)->alfa == NULL){
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





   Read_M (conf_PathR, R, conf_NModels);

   printf ("The model has been read \n");


   //here the model is read

   //calculamos los valores del medelo de las tres ventanas
   //esto hay que meterlo dentro de la estructira de punteros

   Ncoefficients= (pModel)->coef_vector; //Easy code reading. Model is regular so first pModel is used

   ModelParts2= malloc(sizeof(double * )* Ncoefficients); //Allocate space with coefficients size and # models
   if (ModelParts2 == NULL){
      free(pSignal);
      printf ("Not space for ModelParts data array \n");
      exit(0);
   }
   for (i=0;i<Ncoefficients;i++){
	   *(ModelParts2+i)= malloc(sizeof(double ) * conf_NModels);
	   if (*(ModelParts2+i) == NULL){
	      free(pSignal);
	      printf ("Not space for ModelParts data array \n");
	      exit(0);
	   }
   }

   for(j=0;j<Nsignals;j++){

	   for (i=0;i<conf_NModels;i++){ //mean values
	        *(*(ModelParts2+(j*2))+i)= Mean (*((pSignal+j)->pM+i), conf_Npoints);
	   }


	   for (i=0;i<conf_NModels;i++){ //FFT desv.
	        pBuffer= *((pSignal+j)->pM+i);
	        for (t=0;t<conf_Npoints;t++){
	          infft[t][0]=*(pBuffer+t); infft[t][1]=0;
	        }

	        fft(conf_Npoints, infft, outfft);

	        Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

			*(*(ModelParts2+((j*2)+1))+i)= Desv (resabs, conf_Npoints/2);

	      }



   }


  D= (double *) malloc(sizeof(double)* conf_NModels); //Allocate space for distance vectors
  ColumnModel= (double *) malloc(sizeof(double)* Ncoefficients); //Allocate space for auxiliary distance vectors


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

  for(j=0;j<conf_NModels;j++){
	  for (i=0;i<Ncoefficients;i++){
	    *(ColumnModel+i)= *(*(ModelParts2+i)+j);
	  }
	  *(D+j)= distance (ColumnModel, pModel+j);
  }




  //D[2]=  distance (ColumnModel, pModel+2);
  //  printf ("Distancia: %.10f \n",D[0]);
  //  printf ("Distancia: %.10f \n",D[1]);
  // printf ("Distancia: %.10f \n",D[2]);

  //  printf("tR: %f \t %f  \t  %f  \t %f \n",*((pSignal)->pTimeR+31),*((pSignal+2)->pTimeR+31),*((pSignal+3)->pTimeR+31),*((pSignal+7)->pTimeR+31));
  //fprintf(filelog,"tiempo \t S1 \t S2 \t S3 \t S4 \t S5 \t S6 \t S7 \t S8 \t S9 \t S10 \t S11 \t S12 \t S13 \t S14 \t D3 \t D2 \t D1 \t R\n");
  //fprintf(filelog,"tiempo \t S1 \t\t S2 \t S3 \t\t S4 \t\t S5 \t\t S6 \t S7 \t\t S8 \tS9 \t\t S10 \t S11 \t S12 \t S13 \t S14 \t D3 \t   D2 \t    D1 \t      R\n");
  fprintf(filelog,"tiempo \t S1          S2          S3          S4         S5          S6         S7         S8         S9         S10       S11       S12        S13      S14         D3         D2         D1        R\n");

  // printf ("\n Reading distance \n");
  //  fprintf(filelog,"vector D: %.10f \t %.10f \t %.10f \n",D[2],D[1],D[0]);

  Result_Model= (double) ((D[2]*R[2] + D[1]*R[3] + D[0]*R[4] +R[1])/R[0]);

  Result_Model= 0.0;

  for(j=0;j<conf_NModels;j++){
	  Result_Model += (double) (D[j] * R[conf_NModels+1-j]);
  }


  Result_Model= (double)((Result_Model + R[1])/R[0]);
  //Result_Model= (double) ((D[2]*R[2] + D[1]*R[3] + D[0]*R[4] +R[1])/R[0]);
  //printf("Resultado D2: %.10f  D1: %.10f  D0: %.10f \n",D[2],D[1],D[0]);
  // fprintf(filelog,"%.10f \n",Result_Model);

  //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][2],ModelParts[1][2],ModelParts[2][2],ModelParts[3][2],ModelParts[4][2],ModelParts[5][2],ModelParts[6][2],ModelParts[7][2],ModelParts[8][2],ModelParts[9][2],ModelParts[10][2],D[2],D[1],D[0],Result_Model);
  //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][1],ModelParts[1][1],ModelParts[2][1],ModelParts[3][1],ModelParts[4][1],ModelParts[5][1],ModelParts[6][1],ModelParts[7][1],ModelParts[8][1],ModelParts[9][1],ModelParts[10][1],D[2],D[1],D[0],Result_Model);
  //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f  %.7f %.7f %.7f  %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),*(*(ModelParts2+0)+0),*(*(ModelParts2+1)+0),*(*(ModelParts2+2)+0),*(*(ModelParts2+3)+0),*(*(ModelParts2+4)+0),*(*(ModelParts2+5)+0),*(*(ModelParts2+6)+0),*(*(ModelParts2+7)+0),*(*(ModelParts2+8)+0),*(*(ModelParts2+9)+0),*(*(ModelParts2+10)+0),*(*(ModelParts2+11)+0),*(*(ModelParts2+12)+0),*(*(ModelParts2+13)+0),D[2],D[1],D[0],Result_Model);

    fprintf(filelog,"%.3f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f\n",*((pSignal)->pTimeR+31),*(*(ModelParts2+0)+0),*(*(ModelParts2+1)+0),*(*(ModelParts2+2)+0),*(*(ModelParts2+3)+0),*(*(ModelParts2+4)+0),*(*(ModelParts2+5)+0),*(*(ModelParts2+6)+0),*(*(ModelParts2+7)+0),*(*(ModelParts2+8)+0),*(*(ModelParts2+9)+0),*(*(ModelParts2+10)+0),*(*(ModelParts2+11)+0),*(*(ModelParts2+12)+0),*(*(ModelParts2+13)+0),D[2],D[1],D[0],Result_Model);

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
	 ndummy=iwindow1(&iWindow,conf_NModels); //take where to allocate the new buffer window
       j=bufferFree(&iWindow,conf_NModels);

	  for(z=0;z < conf_Npoints;z++){

		  for (i=0;i<Nsignals;i++){ //Resampling for raw signals only =====>  se pasa de conf_Nsignals a Nsignals

			   t=resampling2 ( *(t0+i), conf_Sampling, z, tini,*((pSignal+i)->pM+j) ,(pSignal+i));
		  }
		  tini += conf_Sampling;

	  }

	  for (i=0;i<Nsignals;i++){ //Resampling for raw signals only =====>  se pasa de conf_Nsignals a Nsignals
		  txt[0]= NULL;

		  #ifdef DEBUGLEVEL1
			  sprintf(txt,"Resampleada_%s.txt",(pSignal+i)->name);

			  salvaResampling ( *((pSignal+i)->pM+j), (pSignal+i)->pTimeR, txt, (pSignal+i)->Npoints);
		  #endif

	  }

	  k++;
	  if (k>53)
		  k=105;

     //Here the raw signal are reading and ready to be used

   //calculamos los valores del medelo de las tres ventanas
   //esto hay que meterlo dentro de la estructira de punteros


  for(z=0;z<Nsignals;z++){

     	*(*(ModelParts2+(z*2))+j)= Mean (*((pSignal+z)->pM+j), conf_Npoints);

     	//des ipla
     	pBuffer= *((pSignal+z)->pM+j);
     	for (t=0;t<conf_Npoints;t++){
     		infft[t][0]=*(pBuffer+t); infft[t][1]=0;
     	}

     	fft(conf_Npoints, infft, outfft);

     	Absolute (outfft, resabs, conf_Npoints/2); //Only firts part of FFT

     	//ModelParts[3][i]= Desv (resabs, conf_Npoints/2);
     	*(*(ModelParts2+((z*2)+1))+j)= Desv (resabs, conf_Npoints/2);

   }

  /* for(kk=0;j<conf_NModels;j++){
	   for (i=0;i<Ncoefficients;i++){
		   *(ColumnModel+i)= *(*(ModelParts2+i)+j);
	   }
	   *(D+j)= distance (ColumnModel, pModel+j);
   }*/

   for (i=0;i<Ncoefficients;i++){
	   *(ColumnModel+i)= *(*(ModelParts2+i)+iWindow);
   }


   //  D[0]=  distance (ColumnModel, pModel+iWindow);
   //D[0]=  distance (ColumnModel, pModel); //El modelo siempre es el mismo que D es decir 0 en este caso
   *(D+0)= distance (ColumnModel, pModel);

   ndummy=iwindow_1(&iWindow,conf_NModels);

   for (i=0;i<Ncoefficients;i++){
	   *(ColumnModel+i)= *(*(ModelParts2+i)+ndummy);
   }


  //for (i=0;i<11;i++){
  //  ColumnModel[i]= ModelParts[i][ndummy];
  //}

  // D[2]=  distance (ColumnModel, pModel+ndummy);
  //D[2]=  distance (ColumnModel, pModel+2); //El modelo siempre es el mismo que D es decir 2 en este caso
  *(D+2)= distance (ColumnModel, pModel+2);

  ndummy2=ndummy;

  ndummy= iwindow_1(&ndummy,conf_NModels);

  for (i=0;i<Ncoefficients;i++){
	  *(ColumnModel+i)= *(*(ModelParts2+i)+ndummy);
  }

//  for (i=0;i<11;i++){
//    ColumnModel[i]= ModelParts[i][ndummy];
//  }
  // D[1]=  distance (ColumnModel, pModel+ndummy);
  //D[1]=  distance (ColumnModel, pModel+1); //El modelo siempre es el mismo que D es decir 1 en este caso
  *(D+1)= distance (ColumnModel, pModel+1);




  //  fprintf(filelog,"vector D: %.10f \t %.10f \t %.10f \n",D[2],D[1],D[0]);

   //Result_Model= (double) ((D[2]*R[2] + D[1]*R[3] + D[0]*R[4] +R[1])/R[0]);//ok
   Result_Model= (double) ((*(D+2)*R[2] + *(D+1)*R[3] + *(D+0)*R[4] +R[1])/R[0]);//ok
	#ifdef DEBUGLEVEL2
	   //   Result_Model= D[2]*R[3] + D[1]*R[2] + D[0]*R[1] +R[0];
	   //   printf("Result_Model: %.5f \n",Result_Model);
	   // fprintf(filelog,"%.10f \n",Result_Model);
		#ifdef DEBUGWINDOW
		   //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][ndummy2],ModelParts[1][ndummy2],ModelParts[2][ndummy2],ModelParts[3][ndummy2],ModelParts[4][ndummy2],ModelParts[5][ndummy2],ModelParts[6][ndummy2],ModelParts[7][ndummy2],ModelParts[8][ndummy2],ModelParts[9][ndummy2],ModelParts[10][ndummy2],D[2],D[1],D[0],Result_Model);
		   //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][ndummy],ModelParts[1][ndummy],ModelParts[2][ndummy],ModelParts[3][ndummy],ModelParts[4][ndummy],ModelParts[5][ndummy],ModelParts[6][ndummy],ModelParts[7][ndummy],ModelParts[8][ndummy],ModelParts[9][ndummy],ModelParts[10][ndummy],D[2],D[1],D[0],Result_Model);
		   fprintf(filelog,"%.3f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f\n",*((pSignal)->pTimeR+31),*(*(ModelParts2+0)+ndummy2),*(*(ModelParts2+1)+ndummy2),*(*(ModelParts2+2)+ndummy2),*(*(ModelParts2+3)+ndummy2),*(*(ModelParts2+4)+ndummy2),*(*(ModelParts2+5)+ndummy2),*(*(ModelParts2+6)+ndummy2),*(*(ModelParts2+7)+ndummy2),*(*(ModelParts2+8)+ndummy2),*(*(ModelParts2+9)+ndummy2),*(*(ModelParts2+10)+ndummy2),*(*(ModelParts2+11)+ndummy2),*(*(ModelParts2+12)+ndummy2),*(*(ModelParts2+13)+ndummy2),D[2],D[1],D[0],Result_Model);
		   fprintf(filelog,"%.3f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f\n",*((pSignal)->pTimeR+31),*(*(ModelParts2+0)+ndummy),*(*(ModelParts2+1)+ndummy),*(*(ModelParts2+2)+ndummy),*(*(ModelParts2+3)+ndummy),*(*(ModelParts2+4)+ndummy),*(*(ModelParts2+5)+ndummy),*(*(ModelParts2+6)+ndummy),*(*(ModelParts2+7)+ndummy),*(*(ModelParts2+8)+ndummy),*(*(ModelParts2+9)+ndummy),*(*(ModelParts2+10)+ndummy),*(*(ModelParts2+11)+ndummy),*(*(ModelParts2+12)+ndummy),*(*(ModelParts2+13)+ndummy),D[2],D[1],D[0],Result_Model);
	   #endif
	   fprintf(filelog,"%.3f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f  %.7f\n",*((pSignal)->pTimeR+31),*(*(ModelParts2+0)+iWindow),*(*(ModelParts2+1)+iWindow),*(*(ModelParts2+2)+iWindow),*(*(ModelParts2+3)+iWindow),*(*(ModelParts2+4)+iWindow),*(*(ModelParts2+5)+iWindow),*(*(ModelParts2+6)+iWindow),*(*(ModelParts2+7)+iWindow),*(*(ModelParts2+8)+iWindow),*(*(ModelParts2+9)+iWindow),*(*(ModelParts2+10)+iWindow),*(*(ModelParts2+11)+iWindow),*(*(ModelParts2+12)+iWindow),*(*(ModelParts2+13)+iWindow),D[2],D[1],D[0],Result_Model);
	   //fprintf(filelog,"%.3f\t%.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f   %.7f %.7f %.7f %.7f %.7f\n",*((pSignal)->pTimeR+31),ModelParts[0][iWindow],ModelParts[1][iWindow],ModelParts[2][iWindow],ModelParts[3][iWindow],ModelParts[4][iWindow],ModelParts[5][iWindow],ModelParts[6][iWindow],ModelParts[7][iWindow],ModelParts[8][iWindow],ModelParts[9][iWindow],ModelParts[10][iWindow],D[2],D[1],D[0],Result_Model);

	   //printf ("ventana iwindows: %d  ndummy: %d  ndummy2:  %d  \n", iWindow, ndummy, ndummy2);
	#endif


   if (Result_Model > 0){
     printf("\n ********* Disruption in %d at  t: %f       *********** \n",shotNumber,*((pSignal)->pTimeR+31));
     fprintf (fileRes,"%d \t %s \t %.3f \n",shotNumber,"+1",*((pSignal)->pTimeR+31));
     finish=TRUE;
   }
   else{
//	   if (tini>finaltime || (torigin + 0.010)){
//	   if (tini>dummy){
	   if (tini>finaltime ){
		     printf("\n ********* Signal %d  end NO DISRUPTION       *********** \n",shotNumber ,*((pSignal)->pTimeR+31));
		     fprintf (fileRes,"%d \t %s \t %s \n",shotNumber,"-1","-----");	//If no disruption the alarm time equal 0
		     finish=TRUE;
	   }
   }




   }while(!finish );


	fclose(filelog);
	fclose(fileRes);


    free(pSignal);
    free(pModel);
    free(ColumnModel);




}





void PRINTSIGNAL( signal *entrada){
  int i;

  printf ("Valores comunes Max: %f  Min: %f Pts: %d \n",entrada->Max, entrada->Min, entrada->Npoints);
  printf ("ind \t Tiempo \t\t valor M3 \t valor M2 \t valor M1 \n");
  for (i=0;i<entrada->Npoints;i++){
    printf ("%d \t tr: %f  \t  %f  \t  %f \t  %f \n", i, *(entrada->pTimeR+i),*(entrada->pM+2+i),*(entrada->pM+1+i),*(entrada->pM+i));
  }

}


