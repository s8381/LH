#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <err.h>
#include <string.h>
#include <omp.h>
#include "dot_product.h"
#include "randn.h"
#include "LJ.h"
#include "LJ_shift.h"
//#include "C_remove.h"

#define N 500
#define PI 3.141593
#define k_B 1.380649e-23
#define sig 0.35
#define kappa 1 / 0.35
#define Rs 2
#define Rd 5.5
const float epsilon = 10;
// number of particles
int A_ad_num = 0;
int B_ad_num = 0;
int C_ad_num = 0;
int D_ad_num = 0;
int C_num = 0;
int C_tran;


   
//random number set
   float angle_rand_N[N][2];
 
// Main code
int main(void){    
    //omp_set_num_threads(2); 
    // variables of the simulation
    // scaling: a = 1    D = 1e+9
    const float dt = 1e-5;
    const float evolu = sqrt(2 * dt);
    const float t = 2e2;
    const int Nstep = ceil(t / dt);
    float D_ad = 1;
    float dt_ad = dt / D_ad;
    float evolu_ad = sqrt(2 * dt_ad);
    int loop = 1;
    int sampling_freq = 1000;
    int i, j, k, l, g, h, o, p, q, m, c, cl, co, target, lcxyz, lcyz;
    int lc[3], mcl[3];
    
    lc[0] = 2 * Rd;
    lc[1] = lc[0];
    lc[2] = lc[1];
    lcyz = lc[1] * lc[2];
    lcxyz = lc[0] * lcyz;
    float angle[2];
 
    // lscl[i] holds the atom index to which the i-th atom points
    int buff_l = 100;
    int ls_buff[buff_l];
    int mc[3], rc[3];
    rc[0] = 1;
    rc[1] = 1;
    rc[2] = 1;

    C_tran = 0;
 
    // 1 = A, 2 = B, 3 = Aad, 4 = Bad, 5 = C, 6 = Cad;
    int label[N];
float P_C = 1;//1.4e-6 * 100;
    float r_ad = Rs + 1.1122 * sig;
    float r_escape = Rs + 2.5 * sig;

    // position matrix of the particles
    /*
    float *c_old, *c_new, *r_c;
    c_old = malloc(sizeof(float) * C_num * 3);
    c_new = malloc(sizeof(float) * C_num * 3);
    list = malloc(sizeof(int) * list_size);
    r_c = malloc(sizeof(float) * C_num);
    */
    int *ls = malloc(sizeof(int) * lcxyz * buff_l);    
    int *head = malloc(sizeof(int) * lcxyz);
    float *ra = malloc(sizeof(float) * N * 3);
    float *rb = malloc(sizeof(float) * N * 3);
    float *r = malloc(sizeof(float) * N);
    float r_temp;
    float rd_temp[3];
    

//Brownian dynamics
for(g= 0; g<loop; g++){

    l = 0;    
    char filename1[128], filename2[128], filename3[128];
    snprintf(filename1, sizeof(char) * 128, "traj_%fkT.txt", epsilon);
    FILE *traj;
    traj = fopen(filename1, "w");
    snprintf(filename2, sizeof(char) * 128, "particle_num_%fkT.txt", epsilon);
    FILE *particle_num;
    particle_num = fopen(filename2, "w");
    //snprintf(filename3, sizeof(char) * 64, "traj_c%f.txt", epsilon);
    //FILE *traj_c;
    //traj_c = fopen(filename3, "w");

    #pragma omp parallel for private(i)
    for(i = 0; i < N; i++){
        if(i % 2 == 0)
            label[i] = 1;
	else
            label[i] = 2;
    }
    
    // Random number generator
    srand(time(NULL));

    
// Latice configuration

    int n_temp = 0;
    int mci[3];
    int lci = Rd + 1;
    float lsd = 0.9;
    float ld = 0.5;
    for(mci[0] = 0; mci[0] < (2 * lci + 1); (mci[0])++){
        for(mci[1] = 0; mci[1] < (2 * lci + 1); (mci[1])++){
            for(mci[2] = 0; mci[2] < (2 * lci + 1); (mci[2])++){
                if(n_temp > (N - 1))
                    goto bd;
                for(k = 0; k < 3 ;k++)
                    ra[n_temp * 3 + k] = (-1 * lci + mci[k]) * lsd;
                r[n_temp] = sqrt(dot_product(ra + n_temp * 3, ra + n_temp * 3));
                if(r[n_temp] < (Rd - ld) && r[n_temp] > (Rs + ld)){
                    n_temp++;
                }
                else
                    continue;
            }
        }
    }

    bd:

    
    while(l < Nstep){
	// clear the force and separation matrix
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
	
    // time evolution for A and B
    #pragma omp parallel for private(j, k)
    for(j = 0; j < N; j++){
        float u_vec[3], fr[3];
	float z_strength;
        switch(label[j]){
	    case 1:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 2:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 3:
		for(k = 0; k < 3 ; k++)
		    u_vec[k] = ra[j * 3 + k] / r[j];
		
		for(k = 0; k < 3; k++){
		    fr[k] = evolu * randn(0,1);
		}
		
		z_strength = dot_product(fr, u_vec);
		
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + z_strength * u_vec[k];
                break;
	    case 4:
		for(k = 0; k < 3 ; k++)
		    u_vec[k] = ra[j * 3 + k] / r[j];
		
		for(k = 0; k < 3; k++){
		    fr[k] = evolu * randn(0,1);
		}
		
		z_strength = dot_product(fr, u_vec);
		
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + z_strength * u_vec[k];
                break;
	    case 5:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 6:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 7:
		for(k = 0; k < 3 ; k++)
		    u_vec[k] = ra[j * 3 + k] / r[j];
		
		for(k = 0; k < 3; k++){
		    fr[k] = evolu * randn(0,1);
		}
		
		z_strength = dot_product(fr, u_vec);
		
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + z_strength * u_vec[k];
                break;
	    case 8:
		for(k = 0; k < 3 ; k++)
		    u_vec[k] = ra[j * 3 + k] / r[j];
		
		for(k = 0; k < 3; k++){
		    fr[k] = evolu * randn(0,1);
		}
		
		z_strength = dot_product(fr, u_vec);
		
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + LJ_shift(r[j], epsilon) * ra[j * 3 + k] * dt + z_strength * u_vec[k];
                break;
	}
	
	r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
		
	// Inner B. C. & Outer B. C.
	// reflecting B. C.
        if(r[j] > Rd) {
	    for(k = 0; k < 3; k++)
	        rb[j * 3 + k] = ra[j * 3 + k];
	    // C -> A   D -> B when r = Rd
	    switch(label[j]){
		case 5:
		    label[j] = 1;
		    break;
		case 6:
		    label[j] = 2;
		    break;
	    }
	}

	
	// adsorptopn & desorption of A and B
	switch(label[j]){
            case 1:
                if(r[j] < r_ad){
                    label[j] = 3;
		    A_ad_num++;
		}
                break;
                
            case 2:
                if(r[j] < r_ad){
                    label[j] = 4;
		    B_ad_num++;
		}
                break;
            
            case 3:
                if(r[j] > r_escape){
                    label[j] = 1;
		    A_ad_num--;
		}
                break;
                
            case 4:
                if(r[j] > r_escape){
                    label[j] = 2;
		    B_ad_num--;
		}
                break;
		
	    case 5:
	        if(r[j] < r_ad){
		    label[j] = 7;
		    C_ad_num++;
		}
		break;
	    case 6:
	        if(r[j] < r_ad){
		    label[j] = 8;
		    D_ad_num++;
		}
                break;
		
	    case 7:
	        if(r[j] > r_escape){
		    label[j] = 5;
		    C_ad_num--;
		}
		break;
		
	    case 8:
	        if(r[j] > r_escape){
		    label[j] = 6;
		    D_ad_num--;
		}
		break;
        }
    }
    /*
    // time evolution for C
    if(C_num != 0){
	for(j = 0; j < C_num; j++){
		for(k = 0; k < 3; k++)
		    c_new[j * 3 + k] = c_old[j * 3 + k] + LJ_shift(r_c[j], epsilon) * c_old[j * 3 + k] * dt + evolu * randn(0, 1);
		r_c[j] = sqrt(dot_product(c_new + j * 3, c_new + j * 3));
	    
	    // Reflecting B. C. for C
	    if(r_c[j] < Rs)
		for(k = 0; k < 3; k++)
		    c_new[j * 3 + k] = c_old[j * 3 + k];
	    
	    // record which C sould be removed in the next step
	    if(r_c[j] > Rd){
		C_num--;
		list_size++;
		float *temp_list = realloc(list, sizeof(float) * list_size);
                list = temp_list;
		list[list_size - 1] = j;
	    }
	}
        
	if(list_size != 0){
	    printf("%d C will bw removed\n", list_size);
	    for(j = 0; j < list_size; j++)
		printf("%d ", list[j]);
	    printf("\n");
	    float *temp_new = remove_C(C_num, list_size, list, c_new);
            free(c_new);
            c_new = temp_new;
            list_size = 0;
	    float *temp_r_c = realloc(r_c, sizeof(float) * C_num);
	    r_c = temp_r_c;
	    for(j = 0; j < C_num; j++)
		r_c[j] = sqrt(dot_product(c_new + j * 3, c_new + j * 3));
        }
    }
    */
    

        
	// produce C
	for(i = 0; i < N; i++){
            if(label[i] == 3 || label[i] == 4){
	        int mct[3];
                for(k = 0; k < 3; k++)
		    mct[k] = floor(rb[i * 3 + k] / rc[k]) + Rd;
		
                for(mcl[0] = mct[0]-1; mcl[0] <= mct[0]+1; (mcl[0])++)
                    for(mcl[1] = mct[1]-1; mcl[1] <= mct[1]+1; (mcl[1])++)
                        for(mcl[2] = mct[2]-1; mcl[2] <= mct[2]+1; (mcl[2])++){
                            // Calculate the scalar cell index of the neighbor cells
                            cl = ((mcl[2] + lc[2]) % lc[2]) * lcyz + ((mcl[1] + lc[1]) % lc[1]) * lc[2] + ((mcl[0] + lc[0]) % lc[0]);
                            //scan atom j in cell c1
                            for(p = 0; p < head[cl]; p++){
                                j = ls[cl * buff_l + p];
                                if((label[i] == 3 && label[j] == 4) || (label[i] == 4 && label[j] == 3)){
                                    for(k = 0; k < 3; k++)
                                        rd_temp[k] = rb[i * 3 + k] - rb[j * 3 + k];
                                    r_temp = sqrt(dot_product(rd_temp, rd_temp));
				    if((r_temp < 1.1122 * sig) && ((float) rand() / RAND_MAX <= P_C)){
				        C_tran++;
			                /*
			                C_num++;
			                // expand the matrix of C and assign the position
			                float *temp_pc = realloc(c_new, sizeof(float) * C_num * 3);
			                float *temp_r_c = realloc(r_c, sizeof(float) * C_num);
                                        c_new = temp_pc;
                                        r_c = temp_r_c;
		                        for(k = 0; k < 3; k++)
			                    c_new[(C_num - 1) * 3 + k] = (rb[i][k] + rb[j][k]) / 2;
			                r_c[C_num - 1] = sqrt(dot_product(c_new + (C_num - 1) * 3, c_new + (C_num - 1) * 3));
			                */
					switch(label[i]){
					  case 3:
					    label[i] = 7;
					    label[j] = 8;
					    A_ad_num--;
					    B_ad_num--;
					    C_ad_num++;
					    D_ad_num++
					    break;
					  case 4:
					    label[i] = 8;
					    label[j] = 7;
					    A_ad_num--;
					    B_ad_num--;
					    C_ad_num++;
					    D_ad_num++;
					    break;
					}
				    }
                                }
                            }
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
	
	/*
	fprintf(traj_c, "%d\n\n", C_num);
	    for(j = 0; j < C_num; j++)
		fprintf(traj_c, "%f %f %f %f\n", c_new[j * 3], c_new[j * 3 + 1], c_new[j * 3 + 2], r_c[j]);
	*/
    }

    if(l % 10 == 0)
        fprintf(particle_num, "%d %d %d %d %d\n", A_ad_num, B_ad_num, C_tran, C_ad_num, D_ad_num);
 
    // Relocation for next loop
    /*
    float *temp_old = realloc(c_old, sizeof(float) * C_num * 3);
    c_old = temp_old;
    memcpy(c_old, c_new, C_num * sizeof(float) * 3);
    */
    
    memcpy(ra, rb, N * sizeof(float) * 3);
    
    

    l++;
    }
    /*
    free(c_old);
    free(c_new);
    free(list);
    free(r_c);
    fclose(traj_c);
    */
    fclose(traj);
    fclose(particle_num);
    }
    free(ls);
    free(head);
    free(ra);
    free(rb);
    free(r);
    return 0;
    }
