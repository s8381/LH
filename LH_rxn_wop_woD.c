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
#include "tengential_harmonic.h"


#define N 500
#define PI 3.141593
float sig = 0.35;
float sig_A, sig_B;
float Rs = 2;
float Rd = 5.5;
float epsilon_pp = 0.1;
float epsilon = 10;
float kh = 1000;
int A_ad_num = 0;
int B_ad_num = 0;
int C_num = 0;


// Main code
int main(void){    
    // variables of the simulation
    // scaling: a = 1    D = 1e+9
    const float dt = 1e-5;
    const float evolu = sqrt(2 * dt);
    const float t = 2e3;
    const int Nstep = ceil(t / dt);
    int loop = 1;
    int sampling_freq = 1000;
    int i, j, k, l, g, h, o, p, m, c, cl, lcxyz, lcyz;
    int lc[3], mcl[3];

    sig_A = sig;
    sig_B = sig;    
    lc[0] = 2 * Rd;
    lc[1] = lc[0];
    lc[2] = lc[1];
    lcyz = lc[1] * lc[2];
    lcxyz = lc[0] * lcyz;
    float angle[2];
 
    // lscl[i] holds the atom index to which the i-th atom points  
    int mc[3], rc[3];
    int buff_l = 100;
    int ls_buff[buff_l];
    rc[0] = 1;
    rc[1] = 1;
    rc[2] = 1;
    
    // position matrix of the particles
    int *ls = malloc(sizeof(int) * lcxyz * buff_l);    
    int *head = malloc(sizeof(int) * lcxyz);
    float *ra = malloc(sizeof(float) * N * 3);
    float *rb = malloc(sizeof(float) * N * 3);
    float *r = malloc(sizeof(float) * N);
    float *f = malloc(sizeof(float) * N * 3);
    float *rs = malloc(sizeof(float) * N * 3);
    float r_temp;
    float rd_temp[3];
    float f_buff;
    float fd_buff[3];
 
    // 1 = A, 2 = B, 3 = Aad, 4 = Bad, 5 = C, 6 = Cad;
    int label[N];
    int en = 15;

    // r_inflection = 1.24446
    float P_C = 1;
    float r_ad = Rs + 1.24446 * sig;
    float r_escape = Rs + 1.24446 * sig;
    
// Brownian dynamics
for(g= 0; g < 1; g++){
    l = 0;
    // create traj file
    char filename1[128], filename2[128], filename3[128];
    snprintf(filename1, sizeof(char) * 128, "/media/s8381/HD1/s8381/LH/rxn/r_inflection/wop_woD/k_inf/traj_k_inf_%.1f_long.txt", epsilon);
    FILE *traj;
    traj = fopen(filename1, "w");
    snprintf(filename2, sizeof(char) * 128, "/media/s8381/HD1/s8381/LH/rxn/r_inflection/wop_woD/k_inf/particle_num_k_inf_%.1f_long.txt", epsilon);
    FILE *particle_num;        
    particle_num = fopen(filename2, "w");
    /*
    snprintf(filename3, sizeof(char) * 64, "traj_c%f.txt", epsilon);
    FILE *traj_c;
    traj_c = fopen(filename3, "w");
    */
    
    for(i = 0; i < N / 2; i++)
        label[i] = 1;
    for(i = N / 2; i < N; i++)
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
        // clear the force and header matrix
        memset(f,0,sizeof(float) * N * 3);
        memset(head, 0, sizeof(int) * lcxyz);
        memset(ls, -1, sizeof(int) * lcxyz * buff_l);
        
        // Scan over atoms to construct headers, head and linked list
        for(i = 0; i < N; i++){
            // Calculate the cell vector of atom o
            for(k = 0; k < 3; k++)
                mc[k] = floor(ra[i * 3 + k] / rc[k]) + Rd;
            // translate the cell vector of atom o to cell index
            c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];
            // put particle i into cell c
            ls[c * buff_l + head[c]] = i;
            // increment of head 
            head[c]++;
        }
        
        // Scan inner cells
        for(mc[0] = 0; mc[0] < lc[0]; (mc[0])++){
            for(mc[1] = 0; mc[1] < lc[1]; (mc[1])++){
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
                                            for(k = 0; k < 3; k++)
                                                rd_temp[k] = ra[i * 3 + k] - ra[j * 3 + k];
                                            r_temp = sqrt(dot_product(rd_temp, rd_temp));
                                            float LJ_buff = LJ(r_temp);
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
        float u_vec[3];
        float *harmonic_trapping;
        switch(label[j]){
            case 1:
                for(k = 0; k < 3; k++)
                    rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                    break;
                
            case 2:
                for(k = 0; k < 3; k++)
                    rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                    break;
                
            case 3:
                for(k = 0; k < 3; k++){
                    u_vec[k] = rs[j * 3 + k] / sqrt(dot_product(rs + j * 3, rs + j * 3));
                }
                harmonic_trapping = t_harmonic(u_vec, ra + j * 3);	
            
                for(k = 0; k < 3; k++)
                    rb[j * 3 + k] = ra[j * 3 + k] + (LJ_shift(r[j]) * ra[j * 3 + k] + f[j * 3 + k] + harmonic_trapping[k]) * dt + evolu * randn(0, 1);
                free(harmonic_trapping);
                break;
                
            case 4:
                for(k = 0; k < 3; k++){
                    u_vec[k] = rs[j * 3 + k] / sqrt(dot_product(rs + j * 3, rs + j * 3));
                }
                harmonic_trapping = t_harmonic(u_vec, ra + j * 3);
            
                for(k = 0; k < 3; k++)
                    rb[j * 3 + k] = ra[j * 3 + k] + (LJ_shift(r[j]) * ra[j * 3 + k] + f[j * 3 + k] + harmonic_trapping[k]) * dt + evolu * randn(0, 1);
                free(harmonic_trapping);
                break;
        }
         

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
                for(k = 0; k < 3; k++)
                    rs[j * 3 + k] = rb[j * 3 + k];
                    A_ad_num++;
                }
                break;
                
            case 2:
                if(r[j] < r_ad){
                    label[j] = 4;
                for(k = 0; k < 3; k++)
                    rs[j * 3 + k] = rb[j * 3 + k];
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
        }
    }
    
	// produce C
	for(i = 0; i < N; i++){
	    for(j = 0; j < N; j++){
		if(i == j)
		    continue;
		if((label[i] == 3 && label[j] == 4) || (label[i] == 4 && label[j] == 3)){
		    for(k = 0; k < 3; k ++)
			rd_temp[k] = rb[i * 3 + k] - rb[j * 3 + k];
		    r_temp = sqrt(dot_product(rd_temp, rd_temp)); 
		    
		    if(r_temp < 1.1122 * sig && (float) rand() / RAND_MAX < P_C){
		        // expand the matrix of C and assign the position
		        C_num++;
			// change the labels of AB pair
		        switch(label[i]){
			  case 3:
			    label[i] = 1;
			    label[j] = 2;
			    A_ad_num--;
			    B_ad_num--;
			    break;
			  case 4:
			    label[i] = 2;
			    label[j] = 1;
			    A_ad_num--;
			    B_ad_num--;
			    break;
			}
			
			// let AB pair vanish on the surface and re-insert A and B at the boundary
			int overlapA = 1;
			while(overlapA == 1){
			    reinsertA:;
			    angle[0] = 2 * PI * (float) rand() / RAND_MAX;
			    angle[1] = acos(2 * (float) rand() / RAND_MAX - 1);
			    rb[i * 3] = Rd * sin(angle[0]) * cos(angle[1]);
			    rb[i * 3 + 1] = Rd * sin(angle[0]) * sin(angle[1]);
			    rb[i * 3 + 2] = Rd * cos(angle[0]);
			    
			    // Calculate the cell vector of atom j
			    for(k = 0; k < 3; k++)
				mc[k] = floor(rb[i * 3 + k] / rc[k]) + Rd;
			    // translate the cell vector of atom o to cell index
			    c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];
			    for(mcl[0] = mc[0]-1; mcl[0] <= mc[0]+1; (mcl[0])++)
				for(mcl[1] = mc[1]-1; mcl[1] <= mc[1]+1; (mcl[1])++)
				    for(mcl[2] = mc[2]-1; mcl[2] <= mc[2]+1; (mcl[2])++){
					// Calculate the scalar cell index of the neighbor cells
					cl = ((mcl[2] + lc[2]) % lc[2]) * lcyz + ((mcl[1] + lc[1]) % lc[1]) * lc[2] + ((mcl[0] + lc[0]) % lc[0]);
					//scan atom i in cell c1
					for(p = 0; p < head[cl]; p++){
                                            m = ls[cl * buff_l + p];
                                            if(m == j)
			                        continue;
                                            else{
                                                // check the separation of pair i-j
					        for(k = 0; k < 3; k++)
					            rd_temp[k] = rb[m * 3 + k] - rb[i * 3 + k];
					        r_temp = sqrt(dot_product(rd_temp,rd_temp));
					        if(r_temp < 0.5)
						    goto reinsertA;
                                            }
                                        }
				    }
			    overlapA = 0;
			}
			
			r[i] = sqrt(dot_product(rb + i * 3, rb + i * 3));
                        // Scan over atoms to construct headers, head and linked list
                        for(p = 0; p < N; p++){
                            // Calculate the cell vector of atom o
                            for(k = 0; k < 3; k++)
                                mc[k] = floor(rb[p * 3 + k] / rc[k]) + Rd;
                            // translate the cell vector of atom o to cell index
                            c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];
                            // link to the previous occupant
                            ls[c * buff_l + head[c]] = p;
                            // the last one goes to the header
                            head[c]++;                       
                        }
			// let AB pair vanish on the surface and re-insert A and B at the boundary
			int overlapB = 1;
			while(overlapB == 1){
			    reinsertB:;
			    angle[0] = 2 * PI * (float) rand() / RAND_MAX;
			    angle[1] = acos(2 * (float) rand() / RAND_MAX - 1);
			    rb[j * 3] = Rd * sin(angle[0]) * cos(angle[1]);
			    rb[j * 3 + 1] = Rd * sin(angle[0]) * sin(angle[1]);
			    rb[j * 3 + 2] = Rd * cos(angle[0]);
			    
			    // Calculate the cell vector of atom j
			    for(k = 0; k < 3; k++)
				mc[k] = floor(rb[j * 3 + k] / rc[k]) + Rd;
			    // translate the cell vector of atom o to cell index
			    c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];
			    for(mcl[0] = mc[0]-1; mcl[0] <= mc[0]+1; (mcl[0])++)
				for(mcl[1] = mc[1]-1; mcl[1] <= mc[1]+1; (mcl[1])++)
				    for(mcl[2] = mc[2]-1; mcl[2] <= mc[2]+1; (mcl[2])++){
					// Calculate the scalar cell index of the neighbor cells
					cl = ((mcl[2] + lc[2]) % lc[2]) * lcyz + ((mcl[1] + lc[1]) % lc[1]) * lc[2] + ((mcl[0] + lc[0]) % lc[0]);
					//scan atom p in cell c1
                                        for(p = 0; p < head[cl]; p++){
					    m = ls[cl * buff_l + p];
					// check the separation of pair i-j
					    for(k = 0; k < 3; k++)
						rd_temp[k] = rb[m * 3 + k] - rb[j * 3 + k];
					    r_temp = sqrt(dot_product(rd_temp, rd_temp));
					    if(r_temp < 0.5)
						goto reinsertB;
                                        }
				    }
			    overlapB = 0;
			}
			r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
			
                        // Scan over atoms to construct headers, head and linked list
                        for(p = 0; p < N; p++){
                            // Calculate the cell vector of atom o
                            for(k = 0; k < 3; k++)
                                mc[k] = floor(rb[p * 3 + k] / rc[k]) + Rd;
                            // translate the cell vector of atom o to cell index
                            c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];
                            // link to the previous occupant
                            ls[c * buff_l + head[c]] = p;
                            // the last one goes to the header
                            head[c]++;
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
    }

    if(l % 10 == 0)
        fprintf(particle_num, "%d %d %d\n", A_ad_num, B_ad_num, C_num);
 
    // Relocation for next loop
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
    free(f);
    free(rs);
    
    return 0;
    }
