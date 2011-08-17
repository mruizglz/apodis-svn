#ifndef _TOOLS_H_
#define _TOOLS_H_



typedef  struct {
    char name[CHARMAX];  //Name of signal. Acces path to JET BBDD
    double Max;           //Maximum for normalization
    double Min;           //Minimum for normalization
    int Npoints;         //Number of points for model window
    int nSamples;        //Number of samples of the signal
    double *pData;        //Pointer to raw data
    double *pTime;        //pointer to time for raw data signal
    double *pTimeR;       //pointer to resampling times for nearest time window
    double **pM;         //pointer to array of pointer of model windows
    int Normalize;        //the signal has to be Normalize

  } signal;

typedef struct {
  double gamma;
  double bias;
  int coef_vector;
  int nvectors;
  double **data;
  double *alfa;
} model;



int IndexEvent( double *pBufferin, int n_samples, double Threshold, int type);
double IntLin (double xi, double yi, double xf, double yf, double in);
//int  resampling (int index, double resampling, double t0, double *pDataR, signal *wave, int Normalize);
int  resampling (int index, double resampling, double t0, double *pDataR, signal *wave);
double normalize (double Max, double Min, double data);
int  ReadFloatTxt(char *path, double *buffer);
double Mean (double *pData, int npoints);
int ReadNormalizeTxt(char *path, int *buffer);
int ReadModelTxt(char *path, char **buffer);
void  N_vectors(char * path, int *nvectors, int *ncoefficients);
void  M_values (char * path, model *Model);
void  Read_M (char * path, double *M);
void Absolute ( double (*x)[2], double *result, int npoints);
//void Absolute ( double *rea, double *imj, double *result, int npoints);
double Desv (double *pData, int npoints);
//void distance (double (*input)[3], model *Model, int nModels, double *D);
double distance (double *input, model *Model);
int bufferFree (int *indice, int nModelos);
void salvaOutput( double *input, char *path);
int iwindow_1 (int *indice, int nModelos);
double prod_vect (double *input, double *Coeficientes, double bias, int nCoeficientes);
int iwindow1 (int *indice, int nModelos);




#endif
