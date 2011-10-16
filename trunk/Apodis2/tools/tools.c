
/*

Program implemented by Juan Manuel Lopez, Universidad Politecnica de Madrid
july-Augoust 2011
Implementation of Apodis 2

Tools set to Resamplig and Normalize discharge signals


*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iso646.h>
#include "../tipos.h"
#include "tools.h"





//********************************************************
//   This function search for index when a buffer value
//   is greater or lower
//   than the Threshold, and return the index
//   If type= 0 greater type= 1 lower
//*******************************************************


int IndexEvent( double *pBufferin, int n_samples, double Threshold, int type){

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
	 while (*pBufferin <= Threshold && i < n_samples){
	   i++;
	   pBufferin++;
	 }
	 if (i>0)
		 i--;
       break;

    }

     return i;

}

//********************************************************
//   This function makes a linear interpolation
//
//*******************************************************

double IntLin (double xi, double yi, double xf, double yf, double in){

  double m,b;

   m= (yf-yi)/(xf-xi);

  b= yi-m*xi;


  return (m*in+b);




}


//************************************************
//  This function normalize the signal
//***********************************************


double normalize (double Max, double Min, double data){

  // printf("Normalizacion Max: %f Min: %f ",Max,Min);

  return(double) (data-Min)/(Max-Min);



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


//int  resampling (int index, double resampling, double t0, double *pDataR, signal *wave, int Normalize){
int  resampling (int index, double resampling, double t0, double *pDataR, signal *wave){

  int n; // iterations value
  int i;
  int rindex; //index to complete the resampling buffer
  double sampling; //samplind rate of input signal
  double time; //resamplig time sample
  double paddingvalue; //padding Value, last sample
  double intermediate; //Intermediate Buffer for easy code rading

  //  double t0; //To reuse codefrom previous version
  int Normalize; // To reuse code from previous version

  char txt[255];

  sampling = *(wave->pTime+index+1)-*(wave->pTime+index);

  n= wave->Npoints*resampling / sampling; // number of iterations to resampling

  rindex=0;

  // t0= *(wave->pTime+index);
  Normalize= (wave->Normalize);
  //printf ("Resampling  %s  Index: %d Dato: %f Normaliza: %d  \n",wave->name,index,*(wave->pData+index),Normalize);
if (index+n <= wave->nSamples){
  if (n < wave->Npoints){  //

    for (i=0; i< n; i++){
      do{
	time= t0+(((double)rindex+1)*resampling);
	*(wave->pTimeR+rindex)=time;
	intermediate = IntLin (*(wave->pTime+index+i), *(wave->pData+index+i), *(wave->pTime+index+i+1),  *(wave->pData+index+i+1),time);
	if (Normalize){
	  *(pDataR+rindex)= normalize (wave->Max, wave->Min, intermediate);
	  //	printf(" Res: %f \n",*(pDataR+rindex));
	}

	else{
	  *(pDataR+rindex)= intermediate;
//	  printf(" Intermediate: %f \n",*(pDataR+rindex));
        }
        rindex++;
	time= t0+(((double)rindex+1)*resampling); //use one ahead to see the future jump over the limit
      }while (time < *(wave->pTime+index+i+1));
      paddingvalue= *(pDataR+rindex-1);
    }
    for (i=rindex;i< wave->Npoints;i++){ //padding with the last value
      time= t0+(((double)i+1)*resampling);
       *(pDataR+i)= paddingvalue;
       *(wave->pTimeR+i)= time;
    }
  }
   else{ //For real time, see the possibility to use the value of sample-1

    for (rindex=0;rindex< wave->Npoints;rindex++){
      time= t0+(((double)rindex+1)*resampling);
      i=  IndexEvent( wave->pTime+index,wave->Npoints,time,1);
      if (i>0) i--;
      intermediate = IntLin (*(wave->pTime+index+i), *(wave->pData+index+i), *(wave->pTime+index+i+1),  *(wave->pData+index+i+1),time);
      if (Normalize)
    	  *(pDataR+rindex)= normalize (wave->Max, wave->Min, intermediate);
      else
    	  *(pDataR+rindex)= intermediate;
		  *(wave->pTimeR+rindex)= time;
	  }
   }

     txt[0]= NULL;
	 sprintf(txt,"Original_%s.txt",wave->name);

  salvaResampling ( pDataR, wave->pTimeR, txt, wave->Npoints);
}
else
   n=0;
  return (n);
}


//int resampling2 ( *(t0+i), conf_Sampling,tini,*((pSignal+i)->pM+j) ,(pSignal+i)){
int  resampling2 (int index, double resampling, int indexR, double t0, double *pDataR, signal *wave){


	  int rindex; //index to complete the resampling buffer

	  double time; //resamplig time sample
	  double paddingvalue; //padding Value, last sample
	  double intermediate; //Intermediate Buffer for easy code rading

	  //  double t0; //To reuse codefrom previous version
	  int Normalize; // To reuse code from previous version

	  char txt[255];


	  	  // t0= *(wave->pTime+index);
	  Normalize= (wave->Normalize);
	  //printf ("Resampling  %s  Index: %d Dato: %f Normaliza: %d  \n",wave->name,index,*(wave->pData+index),Normalize);
	  rindex= IndexEvent( (wave->pTime+index), wave->nSamples, t0, 1);
	  rindex= rindex+index;
      intermediate = IntLin (*(wave->pTime+rindex), *(wave->pData+rindex), *(wave->pTime+rindex+1),  *(wave->pData+rindex+1),t0);
      if (Normalize){
    	  *(pDataR+indexR)= normalize (wave->Max, wave->Min, intermediate);
    	  //	printf(" Res: %f \n",*(pDataR+rindex));
      }

      else{
    	  *(pDataR+indexR)= intermediate;
    	  //	  printf(" Intermediate: %f \n",*(pDataR+rindex));
      }

	  *(wave->pTimeR+indexR)= t0;



}



/**********************************************
// This fuction read a txt file and convert
// each line to a double.
// path: Pointer to Path name
// buffer: pointer to result's array
// return the number of lines read
**********************************************/



int ReadFloatTxt(char *path, double *buffer){

	FILE *file;


	char buf[1025], *b, *c;
	int line = 0;
	char *endptr;
	double tmp;


	endptr= NULL;


	file = fopen(path, "r");
	if (!file){
	  printf ("File in ReadFloatTxt Not Found \n");
	  exit(EXIT_FAILURE);
	}
	while (fgets(buf, 1024, file) != NULL)
	{
		b=buf;
		*(buffer+line)= (double)strtod(b,&endptr);
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

  return (result);

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
	int i=0;
	int len;
	char *p;

	p= path;
	i=1;
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
	    for (i=1;i<len;i++){
	      if(isspace(buf[i]) && !isspace(buf[i-1]) && ((buf[i] != '\r') &&  (buf[i] != '\n')  )){
		t++;
	      }
	    }
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
	double dummy;

	//printf("Entra \n");
	file = fopen(path, "r");

	if (!file){
	  printf ("File for Model coeficientes Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	}

	fgets(buf, 1024, file);
	b=buf;
	endptr= NULL;
        //printf("b : %s \n",b);
	Model->gamma= (double)strtod(b,&endptr);
	Model->bias= (double)strtod(endptr,&endptr);
	Model->w=(double)strtod(endptr,&endptr);

	//printf("Gamma: %f   bias: %f w: %f \n",Model->gamma,Model->bias,Model->w);

	for (i=0;i<Model->nvectors;i++){
	  fgets(buf, 1024, file);

	  b=buf;
	  endptr= NULL;
	  t=0;
	  Model->data[i][t]= (double)strtod(b,&endptr);
	  for (t=1;t<Model->coef_vector+1;t++){
	     dummy= (double)strtod(endptr,&endptr);

	     if(t==Model->coef_vector){
	       *(Model->alfa+i)= dummy;
	       //	       printf("valor: %f \n",dummy);
	     }
	     else
	       Model->data[i][t]= dummy;
	     // printf("valor: %f \t",dummy);
	  }

	  line++;


	}


}

void  Read_M (char * path, double *M){


	FILE *file;

	char buf[1025], *b;
	int line = 0;
	int t=0;
	int i;
	int len;
	char *endptr;
	double dummy;

	file = fopen(path, "r");

	if (!file){
	  printf ("File for Model coeficientes Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	}

	fgets(buf, 1024, file);
	b=buf;
	endptr= NULL;
	*(M+4)= (double)strtod(b,&endptr);
	*(M+3)= (double)strtod(endptr,&endptr);
	*(M+2)= (double)strtod(endptr,&endptr);
	*(M+1)= (double)strtod(endptr,&endptr);
	*M= (double)strtod(endptr,&endptr);



}

//void Absolute ( double *rea, double *imj, double *result, int npoints){
void Absolute ( double (*x)[2], double *result, int npoints){

  int i;
  double p1,p2,p3;


  for (i=0;i<npoints;i++){
    *(result+i)= sqrt ( pow (x[i][0],2) + pow (x[i][1],2));
   }


}



double Desv (double *pData, int npoints){ //Ojo se puede ir el puntero y podria retornar infinito


  double  result;
  double  dc;

  int i;

  result=0;
  dc=0;

  for (i=1; i< npoints; i++){ //remove the firt element
    dc+= *(pData+i);
   }
  dc= (double) dc/(npoints-1);
  result=0;
  for (i=1; i< npoints; i++){
    result +=  pow ((*(pData+i) - dc),2) ;
  }
  dc = (double) result / (npoints-2);
  result= sqrt(dc);

  return(result);


}


double distance (double *input, model *Model){

   int i,j,t;
   double *ind_vector;
   double *ind_X;
  double dummy;
  double  norm=0;
   double e=0;
   double  sum=0;
   double u=0;



   ind_vector= *(Model->data);
   ind_X= input;

	sum=0;

	for (i=0;i < Model->nvectors; i++){
                norm=0;
		//	printf("vector: %d \n ", i);
		for (t=0;t < Model->coef_vector;t++){
		   dummy= *(ind_vector+t)-*(ind_X+t);

		   //if((i<2) or (i>Model->nvectors-2))
		     //		     printf("Distancia %.10f in: %.10f  Vector X: %.10f t: %d \n",dummy,*(ind_vector+t),*(ind_X+t),t);
		   //  printf("Vector %.10f \t",*(ind_vector+t));
		   //printf("Input %.10f \t",*(ind_X+t));
		   //printf("Resta %.10f \t",dummy);

		   dummy= pow(dummy,2);
		   //printf("Cuadrado : %.10f\t",dummy);
		   norm=norm+dummy;
		   //printf("norma : %.10f\n",norm);
		}
		ind_vector= *(Model->data+1+i);
		//printf("norma : %.10f\n",norm);
 		//norm= pow(norm,2);
		// printf("x2 : %.10f\n",norm);
		norm=norm*-1.00*Model->gamma;
		//printf("-gamma : %.10f\t",norm);
		e=exp(norm);
		//printf("e : %.10f\t",e);
		u= (*(Model->alfa+i)*e);
		//printf("alfa : %.10f\t",*(Model->alfa+i));
		sum=sum + u;
		//printf("sum : %.10f\n",sum);
	}
	return ((double) (sum+Model->bias)/Model->w);

}


int bufferFree (int *indice, int nModelos){

     if(*indice != 0){
       --*(indice);
     }
     else{
       *(indice)= nModelos-1;
     }

     return (*indice);

}

void salvaOutput( double *input, char *path, int points){

	FILE *file;

	char buf[1025], *b;
	int line = 0;
	int t=0;
	int i;
	int len;
	char *endptr;
	double dummy;

	double *ind_X;



	file = fopen(path, "w");

	if (!file){
	  printf ("File for Output Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	}

	  //	  ind_X= *(input+i);
	  for (t=0;t<points;t++){
	    fprintf(file,"%.10f \t\n",*(input+t));
          }


	fclose(file);



}

int iwindow_1 (int *indice, int nModelos){


  return((*indice - 1) < 0 ? nModelos-1 : *indice -1);


}

int iwindow1 (int *indice, int nModelos){


  return((*indice + 1) > nModelos-1 ? 0 : *indice + 1);


}


double prod_vect (double *input, double *Coeficientes, double bias, int nCoeficientes){

  int i;
  double sum;


  sum=0;
  for (i=0;i<nCoeficientes;i++){
    sum+= *(input+i) * *(Coeficientes+i);
  }

  return (sum+bias);



}

void salvaResampling ( double *input, double *tiempo, char *path, int points){

	FILE *file;

	char buf[1025], *b;
	int line = 0;
	int t=0;
	int i;
	int len;
	char *endptr;
	double dummy;

	double *ind_X;

    b=strchr(path,'/');
    if (b!=NULL)
    	*(b)='_';

	file = fopen(path, "a");

	if (!file){
	  printf ("File for Output Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	}

	  //	  ind_X= *(input+i);
	  for (t=0;t<points;t++){
	    fprintf(file,"%.10f \t %.10f \t\n",*(tiempo+t),*(input+t));
          }


	fclose(file);



}

double  tsincro (char * path, double tini){

	FILE *file;

	char buf[1025], *b;
	int line = 0;
	int t=0;
	int i;
	int len;
	char *endptr;
	double dummy;
	double time, value;

	//printf("Entra \n");
	file = fopen(path, "r");

	if (!file){
	  printf ("Proccess file for time syncro Not Found %s \n", path);
	  exit(EXIT_FAILURE);
	}

	time=0;

	 do {

		fgets(buf, 1024, file);
		b=buf;
		endptr= NULL;
	        //printf("b : %s \n",b);
		time= (double)strtod(b,&endptr);
		value= (double)strtod(endptr,&endptr);


	} while (time < tini);

	 return (time);



}



