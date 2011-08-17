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
    float *pM;           //pointer to array of pointer of model windows
   

  } signal;




int IndexEvent( float *pBufferin, int n_samples, float Threshold, int type);
float IntLin (float xi, float yi, float xf, float yf, float in);
int  resampling (int index, double resampling, float t0, float *pDataR, signal *wave, int Normalize);
float normalize (float Max, float Min, float data);
int  ReadFloatTxt(char *path, float *buffer);
double Mean (double *pData, int npoints);






#endif
