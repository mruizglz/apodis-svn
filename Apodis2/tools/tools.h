#ifndef _TOOLS_H_
#define _TOOLS_H_



typedef  struct {
    char name[CHARMAX];  //Name of signal. Acces path to JET BBDD
    float Max;           //Maximum for normalization
    float Min;           //Minimum for normalization
    int Npoints;         //Number of points for model window
    int nSamples;        //Number of samples of the signal
    float *pData;        //Pointer to raw data
    float *pTime;        //pointer to time for raw data signal
    float *pTimeR;       //pointer to resampling times for nearest time window
    float **pM;         //pointer to array of pointer of model windows
    int Normalize;        //the signal has to be Normalize

  } signal;

typedef struct {
  double gamma;
  double bias;
  int coef_vector;
  int nvectors;
  float **data;
  float *alfa;
} model;



int IndexEvent( float *pBufferin, int n_samples, float Threshold, int type);
float IntLin (float xi, float yi, float xf, float yf, float in);
//int  resampling (int index, double resampling, float t0, float *pDataR, signal *wave, int Normalize);
int  resampling (int index, double resampling, float *pDataR, signal *wave);
float normalize (float Max, float Min, float data);
int  ReadFloatTxt(char *path, float *buffer);
float Mean (float *pData, int npoints);
int ReadNormalizeTxt(char *path, int *buffer);
int ReadModelTxt(char *path, char **buffer);
void  N_vectors(char * path, int *nvectors, int *ncoefficients);
void  M_values (char * path, model *Model);
void  Read_M (char * path, float *M);
void Absolute ( double (*x)[2], double *result, int npoints);
double Desv (double *pData, int npoints);
void distance (float (*input)[3], model *Model, int nModels, float *D);





#endif
