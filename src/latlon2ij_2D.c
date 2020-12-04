#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
using namespace std;

double dist_sphere(double lat1, double lon1, double lat2, double lon2,
                   double R);

float mindist(int npx, int npy, float **lon, float **lat, float tlon, float tlat, int &mdx, int &mdy)
{
    float min_dist;
    float md = 1.e13;
    float R = 6378137; /*Earth's radius im m*/
    int d[8][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};
    int x_next = 0, y_next = 0;

    queue<pair<int, int>> q;
    vector<vector<bool>> visited(npy, vector<bool>(npx, false));
    q.push(make_pair(mdx, mdy));

    while (!q.empty())
    {
        int x = q.front().first;
        int y = q.front().second;
        q.pop();
        visited[y][x] = true;
        min_dist = (float)dist_sphere(lat[y][x], lon[y][x], tlat, tlon, R);
        int count = 0;
        for (int i = 0; i < 8; ++i)
        {
            int x1 = x + d[i][0];
            int y1 = y + d[i][1];
            if (x1 >= 0 && x1 < npx && y1 >= 0 && y1 < npy)
            {
                if (!visited[y1][x1])
                {
                    md = dist_sphere(lat[y1][x1], lon[y1][x1], tlat, tlon, R);
                    if (md <= min_dist)
                    {
                        count++;
                        min_dist = md;
                        x_next = x1;
                        y_next = y1;
                    }
                }
            }
        }
        if (count == 0)
        {
            mdx = x;
            mdy = y;
            return min_dist;
        }
        if (min_dist > 10000)
        {
            q.push(make_pair((x_next + 37) % npx, (y_next + 41) % npy));
            //printf("%d %d\n", (x_next + 37) % npx, (y_next + 41) % npy);
        }
        else
        {
            q.push(make_pair(x_next, y_next));
            //printf("within 10: %d %d %f\n", x_next, y_next, min_dist);
        }
    }
    return min_dist;
}

/*
Input:
------
    NX, NY: int
        Number of grid
    grid_fname: string
        Name of grid FILE
    stat_fname: string
        Name of station ascii file with three columns [site_name, lat, lon]
    output_fname: string
        Name of output ascii file with three columns [site_name, idx, idy]
*/
int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        printf("Inputs: NX, NY, grid_fname, stat_fname, output_fname\nAborting!\n");
        return 0;
    }
    int NX = atoi(argv[1]);
    int NY = atoi(argv[2]);
    char *grid_fname = argv[3];
    char *stat_fname = argv[4];
    char *output_fname = argv[5];

    long int NP, k;
    int mdx = 0, mdy = 0; // Initial trial

    double *buff;
    //float **lat, **lon;
    float *lat[NY], *lon[NY];
    FILE *fid_grid, *fid_in, *fid_out;
    float tlat, tlon, md;
    char cname[15];
    int line, nline = 0, i, j;

    // Read the mesh grid
    NP = (long int)NX * NY;
    buff = (double *)calloc(NP * 3, sizeof(double));
    for (i = 0; i < NY; i++)
    {
        lat[i] = (float *)calloc(NX, sizeof(*(lat[0])));
        lon[i] = (float *)calloc(NX, sizeof(*(lon[0])));
    }
    /* or allocate double pointers
    float **lat = (float**) calloc(NY, sizeof(*lat));
    float **lon = (float**) calloc(NY, sizeof(*lon));
    for (i=0; i < NY; i++) {
        lat[i] = (float*)calloc(NX, sizeof(*(lat[0])));
        lon[i] = (float*)calloc(NX, sizeof(*(lon[0])));
    }
    */

    long int result;
    fprintf(stdout, "Reading mesh %s...\n", grid_fname);
    fflush(stdout);
    fid_grid = fopen(grid_fname, "r");
    result = fread(buff, sizeof(double), NP * 3, fid_grid);
    if (result != NP * 3)
        fprintf(stdout, "Failed to read %ld, \
                result=%ld,3 * NP=%ld\n",
                NP * 3, result, NP);
    fflush(stdout);
    for (int i = 0; i < NY; i++)
    {
        for (int j = 0; j < NX; j++)
        {
            k = (long int)i * NX + j;
            //fprintf(stdout, "k=%ld %f %f\n", k, buff[k * 3], buff[k * 3 + 1]);
            lon[i][j] = (float)buff[k * 3];
            lat[i][j] = (float)buff[k * 3 + 1];
        }
    }
    free(buff);
    fclose(fid_grid);
    fprintf(stdout, "Last grid cell is: %f %f\n", lon[NY - 1][NX - 1], lat[NY - 1][NX - 1]);
    fprintf(stdout, " ok.\n");

    // Main loop
    fid_in = fopen(stat_fname, "r");
    fid_out = fopen(output_fname, "w");
    if (fid_in == NULL)
    {
        fprintf(stdout, "No input station file\n");
        fflush(stdout);
    }
    else
    {
        line = fscanf(fid_in, "%s %f %f\n", cname, &tlon, &tlat);
        while (line != EOF)
        {
            md = mindist(NX, NY, lon, lat, tlon, tlat, mdx, mdy);
            printf("\nProcessing site %d, Dist=%f\n", ++nline, md);
            if (md > 10000)
            {
                fprintf(stdout, "Nearest distance: %f > 10km, probably wrong?\n", md);
            }
            fflush(stdout);
            fprintf(fid_out, "%s %d %d\n", cname, mdx, mdy);
            fflush(fid_out);
            line = fscanf(fid_in, "%s %f %f\n", cname, &tlon, &tlat);
        }
    }
    for (i = 0; i < NY; i++)
    {
        free(lon[i]);
        free(lat[i]);
    }
    fclose(fid_in);
    fclose(fid_out);
    fprintf(stdout, " - finished.\n");
    return (0);
}
