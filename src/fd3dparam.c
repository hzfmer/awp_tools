#include "fd3dparam.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_settings(FD3D_param *P, char *in3dfile)
{
    /* open the two configuration files*/
    FILE *fid = fopen(in3dfile, "r");
    if (fid == NULL)
    {
        fprintf(stdout, "%s\n", in3dfile);
        perror("could not open input file");
        exit(1);
    }

    /* read the parameters from IN3D */
    int status = 2;
    char *stat2;
    char param[15], value[100];
    char newline[200];
    int s;

    stat2 = (char *)calloc(2, sizeof(char));
    while (stat2 != NULL)
    {
        stat2 = fgets(newline, 200, fid);
        /*if (stat2 != NULL)
        fprintf(stdout, "%d: newline=%s\n", strlen(stat2), newline);*/
        if (stat2 != NULL)
        {
            if (strlen(stat2) > 1)
            {
                status = sscanf(newline, "%s %s", value, param);
                if (strcmp(param, "NBGX") == 0)
                    sscanf(value, "%d", &P->nbgx);
                if (strcmp(param, "NBGY") == 0)
                    sscanf(value, "%d", &P->nbgy);
                if (strcmp(param, "NBGZ") == 0)
                    sscanf(value, "%d", &P->nbgz);
                if (strcmp(param, "NEDX") == 0)
                    sscanf(value, "%d", &P->nedx);
                if (strcmp(param, "NEDY") == 0)
                    sscanf(value, "%d", &P->nedy);
                if (strcmp(param, "NEDZ") == 0)
                    sscanf(value, "%d", &P->nedz);
                if (strcmp(param, "NSKPX") == 0)
                    sscanf(value, "%d", &P->nskpx);
                if (strcmp(param, "NSKPY") == 0)
                    sscanf(value, "%d", &P->nskpy);
                if (strcmp(param, "NSKPZ") == 0)
                    sscanf(value, "%d", &P->nskpz);
                if (strcmp(param, "DT") == 0)
                    sscanf(value, "%f", &P->dt);
                if (strcmp(param, "TMAX") == 0)
                    sscanf(value, "%f", &P->tmax);
                if (strcmp(param, "NTISKP") == 0)
                    sscanf(value, "%d", &P->ntiskp);
                if (strcmp(param, "SXRGO") == 0)
                    sscanf(value, "\'%s", P->sxrgo);
                if (strcmp(param, "SYRGO") == 0)
                    sscanf(value, "\'%s", P->syrgo);
                if (strcmp(param, "SZRGO") == 0)
                    sscanf(value, "\'%s", P->szrgo);
                if (strcmp(param, "READ_STEP") == 0)
                    sscanf(value, "%d", &P->readstep);
                if (strcmp(param, "WRITE_STEP") == 0)
                    sscanf(value, "%d", &P->writestep);
                if (strcmp(param, "IVELOCITY") == 0)
                    sscanf(value, "%d", &P->itype);
                if (strcmp(param, "SXRGO") == 0)
                {
                    s = strlen(value);
                    value[s - 1] = '\0';
                    memmove(P->sxrgo, value + 1, s - 1);
                }
                if (strcmp(param, "SYRGO") == 0)
                {
                    s = strlen(value);
                    value[s - 1] = '\0';
                    memmove(P->syrgo, value + 1, s - 1);
                }
                if (strcmp(param, "SZRGO") == 0)
                {
                    s = strlen(value);
                    value[s - 1] = '\0';
                    memmove(P->szrgo, value + 1, s - 1);
                }
            }
        }
    }
    fclose(fid);
}
