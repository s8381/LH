#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <err.h>
#include <string.h>
#include <omp.h>
#include "../interactions/dot_product.h"
#include "../interactions/randn.h"
#include "../interactions/LJ.h"
#include "../interactions/LJ_shift.h"
//#include "C_remove.h"

#define N 2000
#define PI 3.141593
#define k_B 1.380649e-23
#define kappa 1 / 0.35
float epsilon_pp = 0.1;
float epsilon = 3;
float sig = 0.35;
float Rs = 2;
float Rd = 5.5;
const float delta_epsilon = 0;

// number of particles
int A_ad_num = 0;
int B_ad_num = 0;
int C_ad_num = 0;
int D_ad_num = 0;
int C_num = 0;
int *list;
int list_size = 0;
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
    const float t = 1e3;
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
    float P_C = 1;//1.368818e-6;
    float r_ad = Rs + 1.24446 * sig;
    float r_escape = Rs + 1.24446 * sig;
    
    float epsilon_new;
    if(epsilon > delta_epsilon)
        epsilon_new = epsilon - delta_epsilon;
    else
        epsilon_new = 0.1;
    
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
    float *f = malloc(sizeof(float) * N * 3);
    float r_temp;
    float rd_temp[3];
    float f_buff;
    float fd_buff[3];
    

//Brownian dynamics
for(g= 0; g<loop; g++){

    l = 0;    
    char filename1[128], filename2[128], filename3[128];
    snprintf(filename1, sizeof(char) * 128, "/media/s8381/HD1/s8381/LH/rxn/wp/fix_ep_%.0f/N_%d/traj_k_inf_%.1f.txt", epsilon, N / 2, epsilon);
    FILE *traj;
    traj = fopen(filename1, "w");
    snprintf(filename2, sizeof(char) * 128, "/media/s8381/HD1/s8381/LH/rxn/wp/fix_ep_%.0f/N_%d/particle_num_k_inf_%.1f.txt", epsilon, N / 2, epsilon);
    FILE *particle_num;
    particle_num = fopen(filename2, "w");
    //snprintf(filename3, sizeof(char) * 64, "rxn_event_%fkT.txt", epsilon);
    //FILE *rxn_event;
    //rxn_event = fopen(filename3, "w");
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
    float lsd = 0.8;
    float ld = 0.4;
    int overlap = 1;

    float Rd_temp;
    float angle_rand[2];
    angle_rand[0] = 2 * PI * (float) rand() / RAND_MAX;
    angle_rand[1] = acos(2 * (float) rand() / RAND_MAX - 1);
    ra[0] = Rd * sin(angle_rand[0]) * cos(angle_rand[1]);
    ra[1] = Rd * sin(angle_rand[0]) * sin(angle_rand[1]);
    ra[2] = Rd * cos(angle_rand[0]);
    r[0] = sqrt(dot_product(ra, ra));
    
    for(j = 1; j < N; j++){
        overlap = 1;
        while(overlap == 1){
            angle_rand[0] = 2 * PI * (float) rand() / RAND_MAX;
            angle_rand[1] = acos(2 * (float) rand() / RAND_MAX - 1);
            Rd_temp = (Rs + sig) + (Rd - 2 * sig - Rs) * (float) rand() / RAND_MAX;
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
					    float LJ_buff;/*
					    if((label[i] == 1 && label[j] == 5) || (label[i] == 5 && label[j] == 1) || (label[i] == 1 && label[j] == 6) || (label[i] == 6 && label[j] == 1) || (label[i] == 2 && label[j] == 5) || (label[i] == 5 && label[j] == 2) || (label[i] == 2 && label[j] == 6) || (label[i] == 6 && label[j] == 2) || (label[i] == 3 && label[j] == 5)|| (label[i] == 5 && label[j] == 3)|| (label[i] == 3 && label[j] == 6)|| (label[i] == 6 && label[j] == 3)|| (label[i] == 4 && label[j] == 5)|| (label[i] == 5 && label[j] == 4)|| (label[i] == 4 && label[j] == 6)|| (label[i] == 6 && label[j] == 4)){
					        LJ_buff = 0;
						for(k = 0; k < 3; k++)
						    rd_temp[k] = 0;
					    }
					    else{*/
                                                for(k = 0; k < 3; k++)
                                                    rd_temp[k] = ra[i * 3 + k] - ra[j * 3 + k];
                                                r_temp = sqrt(dot_product(rd_temp, rd_temp));
					        LJ_buff = LJ(r_temp);
					    //}
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

	/*
	for(i = 0; i < N; i++){
	    // pair distance calculation
            for (j = i+1; j < N; j++){
                for (k = 0; k < 3; k++){
                    rd[i][j][k] = ra[i][k] - ra[j][k];
		    rd[j][i][k] = -rd[i][j][k];
		}
                r_pp[i][j] = sqrt(dot_product(rd[i][j], rd[i][j]));
		r_pp[j][i] = r_pp[i][j];
            
            // LJ pair force calculation (i,j)
	    // calculate the forces from particle j acting on particle i
	        if(i == j)
		    continue;
		else
		    for(k = 0; k < 3; k++){
			f_pp[i][j][k] = LJ(r_pp[i][j], epsilon_pp) * rd[i][j][k];
			f_pp[j][i][k] = - f_pp[i][j][k];
		    }
	    }

	    // summation of the forces acting on particle i
	    for(j = 0; j < N; j++){
		for(k = 0; k < 3 ; k++)
		    f[i][k] += f_pp[i][j][k];
	    }
        }
	*/
	   
    // time evolution for A and B
    #pragma omp parallel for private(j, k)
    for(j = 0; j < N; j++){
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
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt_ad + LJ_shift(r[j]) * ra[j * 3 + k] * dt_ad + evolu_ad * randn(0, 1);
                break;
	    case 4:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt_ad + LJ_shift(r[j]) * ra[j * 3 + k] * dt_ad + evolu_ad * randn(0, 1);
                break;
	    case 5:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 6:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 7:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	    case 8:
	        for(k = 0; k < 3; k++)
	            rb[j * 3 + k] = ra[j * 3 + k] + f[j * 3 + k] * dt + LJ_shift(r[j]) * ra[j * 3 + k] * dt + evolu * randn(0, 1);
                break;
	}

	r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
		
	// Inner B. C. & Outer B. C.
	// reflecting B. C.
        if(r[j] > Rd) {
            if(label[j] == 5 || label[j] == 6){
                int scan = 0;
                
                // calculate the original c-index
                int mcx[3];
                for(k = 0; k < 3; k++)
                    mcx[k] = floor(ra[j * 3 + k] / rc[k]) + Rd;
                co = mcx[2] * lcyz + mcx[1] * lc[0] + mcx[0];
                for(h = 0; h < buff_l; h++){
                    if(ls[co * buff_l + h] == j){
                        target = h;
                        break;
                    }
                }
		
                // reflecting boundary condition	    
                for(k = 0; k < 3; k++)
                    rb[j * 3 + k] = ra[j * 3 + k];
                r[j] = sqrt(dot_product(rb + j * 3, rb + j * 3));
                
                // remove the particle from the cell
                memcpy(ls_buff, ls + co * buff_l + target + 1, sizeof(int) * (head[co] - (target + 1)));
                memcpy(ls + co * buff_l + target, ls_buff, sizeof(int) * (head[co] - (target + 1)));
                head[co]--;
                memset(ls + co * buff_l + head[co], -1, sizeof(int) * (buff_l - head[co]));
                
                reinsert:;
                if(scan == 1){
                    angle[0] = 2 * PI * (double) rand() / RAND_MAX;
                    angle[1] = acos(2 * (double) rand() / RAND_MAX - 1);
                    rb[j * 3] = Rd * sin(angle[0]) * cos(angle[1]);
                    rb[j * 3 + 1] = Rd * sin(angle[0]) * sin(angle[1]);
                    rb[j * 3 + 2] = Rd * cos(angle[0]);
                    ls[c * buff_l + head[c] - 1] = -1;
                    head[c]--;
                }
                
                
                // calculate the new c-index
                for(k = 0; k < 3; k++)
                    mc[k] = floor(rb[j * 3 + k] / rc[k]) + Rd;
                c = mc[2] * lcyz + mc[1] * lc[0] + mc[0];  
                
                // insert the particle to the new cell
                head[c]++;
                ls[c * buff_l + head[c] - 1] = j;
                
                // check are there A or B particles in the neighboring sites
                // scan the neighboring cells
                for(mcl[0] = mc[0]-1; mcl[0] <= mc[0]+1; (mcl[0])++)
                    for(mcl[1] = mc[1]-1; mcl[1] <= mc[1]+1; (mcl[1])++)
                        for(mcl[2] = mc[2]-1; mcl[2] <= mc[2]+1; (mcl[2])++){
                            
                        // Calculate the scalar cell index of the neighbor cells
                        cl = ((mcl[2] + lc[2]) % lc[2]) * lcyz + ((mcl[1] + lc[1]) % lc[1]) * lc[2] + ((mcl[0] + lc[0]) % lc[0]);
			scan = 1;
                        //scan atom p in cell c1
                        for(p = 0; p < head[cl]; p++){
                            q = ls[cl * buff_l + p];
			    if(q == j)
			        continue;
			    else{
				for(k = 0; k < 3; k++)
				    rd_temp[k] = rb[j * 3 + k] - rb[q * 3 + k];
				r_temp = sqrt(dot_product(rd_temp, rd_temp));
				if(r_temp < sig){
				    goto reinsert;
				}
			    }
                        }
		    }

                // C -> A   D -> B when r = Rd
		relabel:;
                switch(label[j]){
                    case 5:
                        label[j] = 1;
                        break;
                    case 6:
                        label[j] = 2;
                        break;
                }
                
            }
            else{
                for(k = 0; k < 3; k++)
                    rb[j * 3 + k] = ra[j * 3 + k];
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
				    if((r_temp < 1.122462 * sig) && ((float) rand() / RAND_MAX <= P_C)){
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
					    D_ad_num++;
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
					//fprintf(rxn_event, "i = %d, j = %d, t = %d>\n", i, j, l);
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
    //fclose(rxn_event);
    }
    free(ls);
    free(head);
    free(ra);
    free(rb);
    free(r);
    free(f);
    return 0;
    }
