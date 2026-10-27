/*------------------------------------------------------------------------------
* stim300.cc : decode STIM300 imu raw measurement data
*
* version : $Revision:$ $Date:$
* history : 2019/04/14  1.0  new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* read STIM300 IMU raw data--------------------------------------------------
 * args:  char* file  I  STIM300 IMU data file
 *        imu_t *imu  O  imu measurement data
 * return: epochs of imu measurement data
 * --------------------------------------------------------------------------*/
extern int readimucsv(const char *file,imu_t *imu)
{
    imud_t imud={0};
    int week_imu;
    double sec_imu;
    double val[8];
    FILE *handle;

    trace(3,"read imu csv:\n");
    if (!(handle=fopen(file,"rb"))) {
        trace(2,"open imu file fail\n");
        return 0;
    }
    char* headerField = new char[166];
    fread(headerField,sizeof(char),166, handle);
    trace(3,"headerField:%s\n",headerField);
    while (fscanf(handle,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &val[0], &val[1], &val[2], &val[3],
                                                              &val[4], &val[5], &val[6], &val[7]) != EOF) {

        week_imu = (int)val[0];
        sec_imu = val[1];
        imud.time=gpst2time(week_imu,sec_imu);
        imud.gyro[0]=val[2];
        imud.gyro[1]=val[3];
        imud.gyro[2]=val[4];
        imud.accl[0]=val[5];
        imud.accl[1]=val[6];
        imud.accl[2]=val[7];
        addimudata(imu,&imud);
    }
    fclose(handle);
    return imu->n;
}
