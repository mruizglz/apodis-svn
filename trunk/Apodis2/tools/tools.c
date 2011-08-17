/*

Program implemented by Juan Manuel Lopez, Universidad Politecnica de Madrid
july-Augoust 2011
Implementation of Apodis 2

Tools set to Resamplig and Normalize discharge signals


*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../tipos.h"
#include "tools.h"





//********************************************************
//   This function search for index when a buffer value
//   is greater or lower
//   than the Threshold, and return the index
//   If type= 0 greater type= 1 lower
//*******************************************************


int IndexEvent( float *pBufferin, int n_samples, float Threshold, int type){

     int i; // Index

    i=0;

    switch (type){

       case 0:
	 while (*pBufferin > Threshold && i < n_samples){
	   i++;
	   pBufferin++;
	 }
       break;

       case 1:
	 while (*pBufferin < Threshold && i < n_samples){
	   i++;
	   pBufferin++;
	 }
       break;

    }

     return i;

}

//********************************************************
//   This function makes a linear interpolation
//   
//*******************************************************

float IntLin (float xi, float yi, float xf, float yf, float in){

  float m,b;

   m= (yf-yi)/(xf-xi);
 
  b= yi-m*xi;
 

  return (m*in+b);




}


//************************************************
//  This function normalize the signal
//***********************************************


float normalize (float Max, float Min, float data){

  // printf("Normalizacion Max: %f Min: %f ",Max,Min);

  return(float) (data-Min)/(Max-Min);



}


//*************************************************
//  This function makes a resample of a signal
//  index: is initial index for the buffer to use 
//  resampling: the resampling desired 
//  t0: Initial value of the Interpolation buffer start 
//  wave: pointer to struct signal 
//  pDataR: pointer to resampling data buffer
//  normalize: Flag to Normalize, = true (1) Yes
//
//  return the process index value for the next iteration
//*************************************************


//int  resampling (int index, double resampling, float t0, float *pDataR, signal *wave, int Normalize){
int  resampling (int index, double resampling, float *pDataR, signal *wave){

  int n; // iterations value
  int i; 
  int rindex; //index to complete the resampling buffer
  double sampling; //samplind rate of input signal
  float time; //resamplig time sample
  float paddingvalue; //padding Value, last sample
  float intermediate; //Intermediate Buffer for easy code rading

  float t0; //To reuse codefrom previous version
  int Normalize; // To reuse code from previous version

  sampling = *(wave->pTime+1)-*(wave->pTime);

  n= wave->Npoints*resampling / sampling; // number of iterations to resampling

  rindex=0;

  t0= *(wave->pTime+index);
  Normalize= (wave->Normalize);
 
  if (n < wave->Npoints){  //
    for (i=0; i< n; i++){
      do{
	time= t0+((float)rindex*resampling);
	*(wave->pTimeR+rindex)=time;
	intermediate = IntLin (*(wave->pTime+index+i), *(wave->pData+index+i), *(wave->pTime+index+i+1),  *(wave->pData+index+i+1),time);
	if (Normalize){
	  *(pDataR+rindex)= normalize (wave->Max, wave->Min, intermediate);
	  //	printf(" Res: %f \n",*(pDataR+rindex));
	}

	else
	  *(pDataR+rindex)= intermediate;
        rindex++;
	time= t0+((float)rindex*resampling); //use one ahead to see the future jump over the limit
      }while (time < *(wave->pTime+index+i+1));
      paddingvalue= *(pDataR+rindex-1);
    }
    for (i=rindex;i< wave->Npoints;i++){ //padding with the last value
	time= t0+((float)i*resampling); 
       *(pDataR+i)= paddingvalue;
       *(wave->pTimeR+i)= time;
    }
  }
   else{ //For real time, see the possibility to use the value of sample-1
       
    for (rindex=0;rindex< wave->Npoints;rindex++){
      time= t0+((float)rindex*resampling); 
      i=  IndexEvent( wave->pTime+index,wave->Npoints,time,1);
	intermediate = IntLin (*(wave->pTime+index+i), *(wave->pData+index+i), *(wave->pTime+index+i+1),  *(wave->pData+index+i+1),time);
	if (Normalize)
	  *(pDataR+rindex)= normalize (wave->Max, wave->Min, intermediate);
	else
	  *(pDataR+rindex)= intermediate;
       *(wave->pTimeR+rindex)= time;  
    }
   }
  return (n);
}


/**********************************************
// This fuction read a txt file and convert
// each line to a float.
// path: Pointer to Path name
// buffer: pointer to result's array
// return the number of lines read
**********************************************/



int ReadFloatTxt(char *path, float *buffer){

	FILE *file;


	char buf[1025], *b, *c;
	int line = 0;
	char *endptr;
	float tmp;


	endptr= NULL;


	file = fopen(path, "r");
	if (!file){
	  printf ("File in ReadFloatTxt Not Found \n");
	  exit(EXIT_FAILURE);
	} 
	while (fgets(buf, 1024, file) != NULL) 
	{
		b=buf;
		*(buffer+line)= (float)strtod(b,&endptr);
		++line;
	
	}
	return (line);

}


/**********************************************************
//   This function calculate the Mean of Buffer
//   *pData: pointer to buffer
//   npoints; number of points in the buffer
**********************************************************/



double Mean (double *pData, int npoints){ //Ojo se puede ir el puntero y podria retornar infinito


  double result;
  int i;

  result=0;

  for (i=0; i< npoints; i++){
    result+= *(pData+i);
  }

  result= (double) result/npoints;


}

/**********************************************
// This fuction read the Normalize txt file 
// and fill the input array with the normalization Flag
// path: Pointer to Path name
// buffer: pointer to result's array
// return the number of lines read
**********************************************/



int ReadNormalizeTxt(char *path, int *buffer){

	FILE *file;


	char buf[1025], *b, *c;
	int line = 0;
	char *endptr;

	endptr= NULL;


	file = fopen(path, "r");
	if (!file){
	  printf ("File for Normalize Not Found \n");
	  exit(EXIT_FAILURE);
	} 
	while (fgets(buf, 1024, file) != NULL) 
	{
		b=buf;
		*(buffer+line)= strtol(b, &endptr, 0); /* Automatic base 10/16/8 switching */
		++line;
	
	}
	return (line);

}


/**********************************************
// This fuction read the Descripcion Model txt file 
// The file has the path to Model vector Files
// path: Pointer to Path name
// buffer: pointer to result's array
// return the number of lines read
**********************************************/



int ReadModelTxt(char *path, char **buffer){

	FILE *file;


	char buf[1025];
	int line = 0;

	file = fopen(path, "r");
	if (!file){
	  printf ("File for Model description Not Found \n");
	  exit(EXIT_FAILURE);
	} 
	while (fgets(buf, 1024, file) != NULL) 
	{
	
	   strcpy(*(buffer+line), buf); 
	   ++line;
	
	}
	return (line);

}

/***********************************************************
// This fuction return de number of coefficients in a vector
// and the number of vectors
************************************************************/


void  N_vectors(char * path, int *nvectors, int *ncoefficients){


	FILE *file;


	char buf[1025];
	int line = 0;
	int t=0;
	int i;
	int len;
	char *p;

	p= path;
	i=0;
	while (*(p+i++)!= '\n');
	*(p+i-1)= '\0';
	
	file = fopen(path, "r");
	if (!file){
	  printf ("File for coeficientes Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	} 
	while (fgets(buf, 1024, file) != NULL) 
	{
	  if(line == 1){
	    len=strlen(buf);
	    for (i=1;i<len;i++)
	      if(isspace(buf[i]) && !isspace(buf[i-1]) && (buf[i] != '\n'))
		t++;
	    }
	  //else do nothing
	   if (buf[0] != '\n')
	     ++line;
	}

	*nvectors= line-1;
	*ncoefficients= t;

	
}


void  M_values (char * path, model *Model){


	FILE *file;

	char buf[1025], *b;
	int line = 0;
	int t=0;
	int i;
	int len;
	char *endptr;
	float dummy;

	file = fopen(path, "r");

	if (!file){
	  printf ("File for Model coeficientes Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	} 

	fgets(buf, 1024, file);
	b=buf;
	endptr= NULL;
	Model->gamma= (float)strtod(b,&endptr);
	Model->bias= (float)strtod(endptr,&endptr);
	printf("Gamma: %f   bias: %f \n",Model->gamma,Model->bias);
	
	for (i=0;i<Model->nvectors;i++){
	  fgets(buf, 1024, file);
 
	  b=buf;
	  endptr= NULL;
	  t=0;
	  Model->data[i][t]= (float)strtod(b,&endptr);
	  for (t=0;t<Model->coef_vector;t++){
	     dummy= (float)strtod(endptr,&endptr);

	     if(t==Model->coef_vector-1){
	       *(Model->alfa+i)= dummy;
	       printf("valor: %f \n",dummy);
	     }
	     else
	       Model->data[i][t]= dummy;
	  }

	  line++;
 
  
	}

	
}
