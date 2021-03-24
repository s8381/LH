#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <err.h>
#include <string.h>
#include "dot_product.h"
#include "randn.h"
#include "LJ.h"
#include "LJ_shift.h"
#include "LJ6.h"

#define N 500
#define PI 3.141593

float Rs = 2;
float Rd = 5.5;
float epsilon_pp = 0.1;
float epsilon = 10;
float sig = 0.35;
// number of particles
int A_ad_num = 0;
int B_ad_num = 0;
int A_ad_tran, B_ad_tran, A_tran, B_tran;

 
// Main code
int main(void){    
    // variables of the simulation
    // scaling: a = 1    D = 1e+9
    const float dt = 1e-5;
    const float evolu = sqrt(2 * dt);
    const float t = 1e3;
    const int Nstep = ceil(t / dt);
    int loop = 1;
    int sampling_freq = 1000;
    int i, j, k, l, g, h, o, p, c, cl, lcxyz, lcyz, buff_l;
    int lc[3], mcl[3], mc[3];
    
    lc[0] = 2 * Rd;
    lc[1] = lc[0];
    lc[2] = lc[1];
    lcyz = lc[1] * lc[2];
    lcxyz = lc[0] * lcyz;
    buff_l = 200;
    
// position matrix of the particles
    int *ls = malloc(sizeof(int) * lcxyz * buff_l);    
    int *head = malloc(sizeof(int) * lcxyz);
    float *ra = malloc(sizeof(float) * N * 3);
    float *rb = malloc(sizeof(float) * N * 3);
    float *r = malloc(sizeof(float) * N);
    float *f = malloc(sizeof(float) * N * 3);
    float r_temp;
    float rd_temp[3];
    float f_buff;
    float fd_buff[3];
    
    // lscl[i] holds the atom index to which the i-th atom points  
    int rc[3];
    rc[0] = 1;
    rc[1] = 1;
    rc[2] = 1;
    
    float rshitf[3], Region[3];
    Region[0] = 2 * Rd;
    Region[1] = Region[0];
    Region[2] = Region[1];   


    A_ad_tran = 0;
    B_ad_tran = 0;
    A_tran = 0;
    B_tran = 0;
 
    // 1 = A, 2 = B, 3 = Aad, 4 = Bad, 5 = C, 6 = Cad;
    int label[N];

    float P_C = 1;
    float r_ad = Rs + 1.1122 * sig;
    float r_escape = Rs + 2.5 * sig;
    
    
// Brownian dynamics
for(g= 0; g<loop; g++){
    l = 0;
    // create traj file
    char filename1[64], filename2[64], filename3[64];
    snprintf(filename1, sizeof(char) * 64, "traj_%fkT_ad.txt", epsilon);
    FILE *traj;
    traj = fopen(filename1, "w");
    snprintf(filename2, sizeof(char) * 64, "particle_num_%fkT_ad.txt", epsilon);
    FILE *particle_num;        
    particle_num = fopen(filename2, "w");

    
    for(i = 0; i < N/2; i++)
        label[i] = 1;
    for(i = N/2; i < N; i++)
        label[i] = 2;
    
    // Random number generator
    srand(time(NULL));


    // Randomly distribute the particles
    float angle_rand[2];
    angle_rand[0] = 2 * PI * (float) rand() / RAND_MAX;
    angle_rand[1] = acos(2 * (float) rand() / RAND_MAX - 1);
    ra[0] = Rd * sin(angle_rand[0]) * cos(angle_rand[1]);
    ra[1] = Rd * sin(angle_rand[0]) * sin(angle_rand[1]);
    ra[2] = Rd * cos(angle_rand[0]);
    r[0] = sqrt(dot_product(ra, ra));
    int overlap = 1;

    float Rd_temp;
    for(j = 1; j < N; j++){
        overlap = 1;
        while(overlap == 1){
            angle_rand[0] = 2 * PI * (float) rand() / RAND_MAX;
            angle_rand[1] = acos(2 * (float) rand() / RAND_MAX - 1);
            Rd_temp = (Rs + sig) + (Rs - 2 * sig) * (float) rand() / RAND_MAX;
            ra[j * 3] = Rd_temp * sin(angle_rand[0]) * cos(angle_rand[1]);
            ra[j * 3 + 1] = Rd_temp * sin(angle_rand[0]) * sin(angle_rand[1]);
            ra[j * 3 + 2] = Rd_temp * cos(angle_rand[0]);
            for(p = 0; p < j; p++){
                for(k = 0; k < 3; k++)
                    rd_temp[k] = ra[j * 3 + k] - ra[p * 3 + k];
                r_temp = sqrt(dot_product(rd_temp, rd_temp));
                if(r_temp < 0.8 * sig){
                    overlap = 1;
                    break;
                }
                else
                    overlap = 0;
            }
        }
        r[j] = sqrt(dot_product(ra + j * 3, ra + j * 3));
    }
    
    
   
    bd:
    while(l < Nstep){
	// clear the force and separation matrix
	memset(f,0,sizeof(float) * N * 3);
	memset(head, 0, sizeof(int) * lcxyz);
	memset(ls, -1, sizeof(int) * lcxyz * buff_l);

	    
	// Scan over atoms to construct headers, head and linked list
	#pragma omp parallel for private(i, k, c)
	for(i = 0; i < N; i++){
	    int mct[3];
	    // Calculate the cell vector of atom i
	    for(k = 0; k < 3; k++)
		mct[k] = floor(ra[i * 3 + k] / rc[k]) + Rd;
	    // translate the cell vector of atom o to cell index
	    c = mct[2] * lcyz + mct[1] * lc[0] + mct[0];
	    // put particle i into cell c
	    ls[c * buff_l + head[c]] = i;
	    // increment of head 
	    head[c]++;
	}

        
        // Scan inner cells
        #pragma parallel for private(mc, c, mcl, cl, k, i, j, o, p, rd_temp, r_temp, fd_buff)
        for(mc[0] = 0; mc[0] < lc[0]; (mc[0])++){
            //#pragma omp parallel for private(mc)
            for(mc[1] = 0; mc[1] < lc[1]; (mc[1])++){
                //#pragma omp parallel for private(mc)
                for(mc[2] = 0; mc[2] < lc[2]; (mc[2])++){
                    // Calculate a cell index
                    c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];
                    // scan the neighbor cells (and itself)
                    for(mcl[0] = mc[0]-1; mcl[0] <= mc[0]+1; (mcl[0])++)
                        for(mcl[1] = mc[1]-1; mcl[1] <= mc[1]+1; (mcl[1])++)
                            for(mcl[2] = mc[2]-1; mcl[2] <= mc[2]+1; (mcl[2])++){
                                // Calculate the scalar cell index of the neighbor cells
                                cl = ((mcl[2] + lc[2]) % lc[2]) * lcyz + ((mcl[1] + lc[1]) % lc[1]) * lc[2] + ((mcl[0] + lc[0]) % lc[0]);
                                // scan atom o in cell c
                                for(o = 0; o < head[c]; o++){
                                    i = ls[c * buff_l + o];
                                    //scan atom p in cell c1
                                    for(p = 0; p < head[cl]; p++){
                                        j = ls[cl * buff_l + p];
                                        if(i < j){
					    float LJ_buff;
                                            for(k = 0; k < 3; k++)
                                            rd_temp[k] = ra[i * 3 + k] - ra[j * 3 + k];
                                            r_temp = sqrt(dot_product(rd_temp, rd_temp));
                                            LJ_buff = LJ(r_temp);
                                            
					    for(k = 0; k < 3; k++){
						fd_buff[k] = LJ_buff * rd_temp[k];
						f[i * 3 + k] += fd_buff[k];
						f[j * 3 + k] -= fd_buff[k];
					    }
                                        }
                                    }
                                }
                            }
                }
	    }
	}

	    
    // time evolution for A and B
    for(j = 0; j < N; j++){
        // vibration of adsorbed aprticles
	for(k = 0; k < 3; k++)
	    rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
	r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
		
	// Inner B. C. & Outer B. C.
	// reflecting B. C.
        if(r[j] < Rs || r[j] > Rd) {
	    for(k = 0; k < 3; k++)
		rb[j * 3 + k] = ra[j * 3 + k];
	    r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
	}

	
	// adsorptopn & desorption of A and B
	switch(label[j]){
            case 1:
                if(r[j] < r_ad){
                    label[j] = 3;
		    A_ad_num++;
		    A_ad_tran++;
		}
                break;
                
            case 2:
                if(r[j] < r_ad){
                    label[j] = 4;
		    B_ad_num++;
		    B_ad_tran++;
		}
                break;
            
            case 3:
                if(r[j] > r_escape){
                    label[j] = 1;
		    A_ad_num--;
		    A_tran++;
		}
                break;
                
            case 4:
                if(r[j] > r_escape){
                    label[j] = 2;
		    B_ad_num--;
		    B_tran++;
		}
                break;
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

    if(l % 10 == 0)
        fprintf(particle_num, "%d %d\n", A_ad_num, B_ad_num);
    

 
    // Relocation for next loop  
    memcpy(ra, rb, N * sizeof(float) * 3);
    
    

    l++;
    }
    fclose(traj);
    fclose(particle_num);
    }

    return 0;
    }
