#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <err.h>
#include <string.h>
#include <omp.h> 
#include "dot_product.h"
#include "randn.h"
#include "LJ_shift.h"

#define N 200
#define PI 3.141593

// global variables
float epsilon = 0.1;
float sig_A = 0.4;
float sig_B = 0.4;
float sig = 0.4;

float Rs = 1;
float Rd = 10;

/* make the electric energy * length dimensionless, the 
       
  e^2 / (4 * pi * epsilon_r * epsilon_0) / (kT * 1 nm)
  1.6^2 * 10^(-38) / (4 * pi * 78 * 8.85 * 10^(-12)) / (1.38 * 10^(-23) * 298 * 10^(-9)) = 0.7176239

*/

// Main code
int main(void){    
    // variables of the simulation
    // scaling: a = 1    D = 1e+9
    const float dt = 1e-4;
    const float evolu = sqrt(2 * dt);
    int loop = 10;
    int sampling_freq = 10000;
    int event_sampling = 10;
    int i, j, k, l, o, p, g, c, cl;
    int lc, lcyz, lcxyz;
    
    // position matrix of the particles
    float *ra = malloc(sizeof(float) * N * 3);
    float *rb = malloc(sizeof(float) * N * 3);
    float *r = malloc(sizeof(float) * N);
    float r_temp;
    float rd_temp[3];
    float fr[3];

 
    // 1 = A, 2 = B
    int label[N];

    // reaction probability
    float P_C = 1;
     
    // reaction cut-off
    float r_ad = Rs + sig_A * 1.122462;
    

    
// Brownian dynamics
for(g = 0; g < 100; g++){
    // create traj file
    char filename2[128];
    snprintf(filename2, sizeof(char) * 128, "survival_%d.txt", g);
    FILE *particle_num;        
    particle_num = fopen(filename2, "w");
    
    // open one trajectory file
    char filename1[128];
    snprintf(filename1, sizeof(char) * 128, "traj_%d.txt", g);
    FILE *traj;
    traj = fopen(filename1, "w");
    
    // first-passage time file
    char filename3[128];
    snprintf(filename3, sizeof(char) * 128, "FPT_%d.txt", g);
    FILE *FPT;
    FPT = fopen(filename3, "w");

    // conversion of A
    float conversion = 0;    
    
    // total number of A
    float total_A = N;
    
    // all particles are A initially
    for(i = 0; i < N; i++)
        label[i] = 1;
    
    // Random number generator
    srand(time(NULL));

    memset(r, 0, sizeof(float) * N);
    
    // Latice configuration
    float ld = 1.5;

    while(r[0] < Rs + ld){
        for(k = 0; k < 3; k++){
            do{
                ra[k] = randn(0,1) * (Rd - ld);
            }while((ra[k] > (Rd - ld)) || (ra[k] < (-Rd + ld)));
        }
        r[0] = sqrt(dot_product(ra, ra));
    }

    for(i = 1; i < N; i++){
        int overlap = 1;
        while(overlap == 1 || r[i] < Rs + ld){
            for(k = 0; k < 3; k++){
                do{
                    ra[i * 3 + k] = randn(0,1) * (Rd - ld);
                }while(ra[i * 3 + k] > (Rd - ld) || ra[i * 3 + k] < (-Rd + ld));
            }
            // check the distance between the existing particles and the newly inserted one
            for(j = 0; j < i; j++){
                float rd_temp[3], r_temp;
                for(k = 0; k < 3; k++)
                    rd_temp[k] = ra[i * 3 + k] - ra[j * 3 + k];
                r_temp = sqrt(dot_product(rd_temp, rd_temp));
                if(r_temp < ld){
                    overlap = 1;
                    break;
                }
                else
                    overlap = 0;
                
           }
           r[i] = sqrt(dot_product(ra + i * 3, ra + i * 3));
        }
    }

    l = 0;
    // stop the simulation when the conversion of HCF is larger than 80%
    while(conversion < 0.8){

    // time evolution for A and B
    #pragma parallel for private(j, k)
    for(j = 0; j < N; j++){
        for(k = 0; k < 3; k++)
            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
        r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));

	// adsorbing B. C. at the sink, A -> B
        if(r[j] <= r_ad && label[j] == 1){
            label[j] = 2;
            total_A--;
            fprintf(FPT, "%d\n", l);
	}

        // periodic B. C. at the edge of the box
        for(k = 0; k < 3; k++){
	    if(rb[j * 3 + k] > Rd){
                rb[j * 3 + k] -= 2 * Rd;
                r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
            }
            else if(rb[j * 3 + k] < -Rd){
                rb[j * 3 + k] += 2 * Rd;
                r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
            }
        }
    }
    
    // Write trajectories
    if(l % sampling_freq == 0) {
        for (j = 0; j < N; j++)
            r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
        fprintf(traj, "%d\n\n", N);
	    for(j = 0; j < N; j++)
		fprintf(traj, "type%d %f %f %f id=%d> \n", label[j], rb[j * 3], rb[j * 3 + 1], rb[j * 3 + 2], j);
    }

    // print the number of reaction events
    if(l % event_sampling == 0)
        fprintf(particle_num, "%0.f\n", total_A);

    
    // rellocate the trajectory matrix
    memcpy(ra, rb, N * sizeof(float) * 3);
    
    // calculate the conversion of A
    conversion = (N - total_A) / N;
    
    l++;
    }
    fclose(particle_num);
    fclose(traj);
    fclose(FPT);
    }
    free(ra);
    free(rb);
    free(r);
}
