#ifndef _CONFIG_H_
#define _CONFIG_H_

extern char *conf_file_name;
extern double   conf_Threshold;
extern double   conf_Sampling;

extern int      conf_Npoints;
extern int      conf_Nsignals;
extern int      conf_NcalculateSignals;
extern int      conf_NModels;




extern char     *conf_PathMax;
extern char     *conf_PathMin;

extern char     *conf_PathD1;
extern char     *conf_PathD2;
extern char     *conf_PathD3;

extern char     *conf_PathR;


/*
extern int   conf_TOutRs485UCms;
extern int   conf_TOutRs485MODms;
extern int   conf_TOutResetAlarmams;
extern char *conf_devRs485UC;
extern char *conf_devRs485MOD;
extern char *conf_devGateWay;
extern int   conf_noDaemon;

extern char *conf_smtpServer;
*/




void read_config_files(void);


#endif
