#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define N 750
#define PI 3.14159265

float dot_product(float *a, float *b){
    float result;
    int p;
    result = 0;
    for(p = 0; p < 3; p++){
        result += *(a + p) * *(b + p);
    }
    
    return result;
}

float pair_distance(float *a, float *b){
    int i;
    float rd[3], pd;
    
    for(i = 0; i < 3; i++)
        rd[i] = *(a + i) - *(b + i);
    pd = sqrt(dot_product(rd, rd));
    
    return pd;
}

float arclength(float *a, float *b){
    float rd[3], al, theta, R_avg;
    
    theta = acos(dot_product(a, b) / sqrt(dot_product(a, a)) / sqrt(dot_product(b, b)));
    R_avg = (sqrt(dot_product(a, a)) + sqrt(dot_product(b, b))) / 2;
    al = R_avg * theta;
    return al;
}

float pd_normalization(float freq, int r){
    float gr_pd;
    float interval = 0.05;
    gr_pd = freq / (4 * PI * r * r * interval * interval * interval);
    return gr_pd;
}

float al_normalization(float freq, int r){
    float gr_al, x, theta, R_avg, interval;
    R_avg = 2 + 1.122462 * 0.35;
    interval = 0.05;
    theta = r / R_avg * interval;
    gr_al = freq / (2 * PI * R_avg * fabs(sin(theta)) * interval);
    return gr_al;
}


/*
Change

1. file names
2. length_of_gr
3. length calculation
4. histogram

*/ 

int main(){

    int i, j, m, n, id, type, frame, frame_AA, frame_BB, frame_AB, frame_CC, frame_DD, frame_CD, ep, en;
    en = 19;
    ep = 10;
    size_t size = 64;
    char s1[64], s2[64], s3[64], s4[64], s5[64], s6[64], s7[64], s8[64], s9[64], s10[64], s11[64], s12[64], s13[64];
    char *b = s1;
    float r_temp[3], ra[N][3], r, rd[3];
    snprintf(s1, sizeof(char) * 64, "../traj_k_inf_%d.0_%d.txt", ep, en);
    snprintf(s2, sizeof(char) * 64, "g33_%dkT_al_%d.txt", ep, en);
    snprintf(s3, sizeof(char) * 64, "g44_%dkT_al_%d.txt", ep, en);
    snprintf(s4, sizeof(char) * 64, "g34_%dkT_al_%d.txt", ep, en);
    snprintf(s5, sizeof(char) * 64, "freq33_%dkT_al_%d.txt", ep, en);
    snprintf(s6, sizeof(char) * 64, "freq44_%dkT_al_%d.txt", ep, en);
    snprintf(s7, sizeof(char) * 64, "freq34_%dkT_al_%d.txt", ep, en);
    
    FILE *traj = fopen(s1, "r");    
    FILE *g33_out = fopen(s2, "w");
    FILE *g44_out = fopen(s3, "w");
    FILE *g34_out = fopen(s4, "w");
    FILE *freq33_out = fopen(s5, "w");
    FILE *freq44_out = fopen(s6, "w");
    FILE *freq34_out = fopen(s7, "w");
    
/*    
    snprintf(s8, sizeof(char) * 64, "g77_%dkT_al.txt", ep);
    snprintf(s9, sizeof(char) * 64, "g88_%dkT_al.txt", ep);
    snprintf(s10, sizeof(char) * 64, "g78_%dkT_al.txt", ep);
    snprintf(s11, sizeof(char) * 64, "freq_77%dkT_al.txt", ep);
    snprintf(s12, sizeof(char) * 64, "freq_88%dkT_al.txt", ep);
    snprintf(s13, sizeof(char) * 64, "freq_78%dkT_al.txt", ep);    

    FILE *g77_out = fopen(s8, "w");
    FILE *g88_out = fopen(s9, "w");
    FILE *g78_out = fopen(s10, "w");
    FILE *freq77_out = fopen(s11, "w");
    FILE *freq88_out = fopen(s12, "w");
    FILE *freq78_out = fopen(s13, "w");
*/
    
    frame = 0;
    frame_AA = 0;
    frame_BB = 0;
    frame_AB = 0;
    
    frame_CC = 0;
    frame_DD = 0;
    frame_CD = 0;
    r = 0;

    float Rs, sigma, L, aL;
    Rs = 2;
    sigma = 0.35;
    L = 2 * Rs + 5 * sigma;
    aL = PI * (Rs + 2.5 * sigma);

    int A_num, B_num;
    A_num = 0;
    B_num = 0;

    // array of the index of adsorbants
    int *A_ad_array, *B_ad_array;
    A_ad_array = malloc(0 * sizeof(int));
    B_ad_array = malloc(0 * sizeof(int));
    

    // array of g(r)
    float interval = 0.05;
    float al_interval = 0.05;
    int  length_of_gr = round(aL / al_interval);
    B_ad_array = malloc(0 * sizeof(int));
    //C_ad_array = malloc(0 * sizeof(int));
    //D_ad_array = malloc(0 * sizeof(int));
    

    // array of g(r)
    float g33[length_of_gr], g44[length_of_gr], g34[length_of_gr];
    float g77[length_of_gr], g88[length_of_gr], g78[length_of_gr];
    float freq33[length_of_gr], freq44[length_of_gr], freq34[length_of_gr];
    float freq77[length_of_gr], freq88[length_of_gr], freq78[length_of_gr];
    float freq_buff33[length_of_gr], freq_buff44[length_of_gr], freq_buff34[length_of_gr];
    float freq_buff77[length_of_gr], freq_buff88[length_of_gr], freq_buff78[length_of_gr];
    int ri;

    // initialization of arrays
    memset(g33, 0, sizeof(float) * length_of_gr);
    memset(g44, 0, sizeof(float) * length_of_gr);
    memset(g34, 0, sizeof(float) * length_of_gr);
    memset(freq33, 0, sizeof(float) * length_of_gr);
    memset(freq44, 0, sizeof(float) * length_of_gr);
    memset(freq34, 0, sizeof(float) * length_of_gr);
    memset(freq_buff33, 0, sizeof(float) * length_of_gr);
    memset(freq_buff44, 0, sizeof(float) * length_of_gr);
    memset(freq_buff34, 0, sizeof(float) * length_of_gr);

    /*
    memset(g77, 0, sizeof(float) * length_of_gr);
    memset(g88, 0, sizeof(float) * length_of_gr);
    memset(g78, 0, sizeof(float) * length_of_gr);
    memset(freq77, 0, sizeof(float) * length_of_gr);
    memset(freq88, 0, sizeof(float) * length_of_gr);
    memset(freq78, 0, sizeof(float) * length_of_gr);
    memset(freq_buff77, 0, sizeof(float) * length_of_gr);
    memset(freq_buff88, 0, sizeof(float) * length_of_gr);
    memset(freq_buff78, 0, sizeof(float) * length_of_gr);
    */

    // read the trajectory file line by line
    while(getline(&b, &size,traj) != -1){
        int read = sscanf(s1, "type%d %f %f %f id=%d>", &type, &r_temp[0], &r_temp[1], &r_temp[2], &id);

        // store the coordinates into matrix ra
        if(read == 5){
                for(i = 0; i < 3; i++)
                    ra[id][i] = r_temp[i];
        
                // count the number of Aad & Bad and store the index of adsorbants
                switch(type){
                    case 3:
                        A_num++;
                        int *A_ad_array_temp = realloc(A_ad_array, sizeof(int) * A_num);
                        A_ad_array = A_ad_array_temp;
                        A_ad_array[A_num - 1] = id;
                        break;

                    case 4:
                        B_num++;
                        int *B_ad_array_temp = realloc(B_ad_array, sizeof(int) * B_num);
                        B_ad_array = B_ad_array_temp;
                        B_ad_array[B_num - 1] = id;
                        break;
    /*
                    case 7:
                        C_num++;
                        int *C_ad_array_temp = realloc(C_ad_array, sizeof(int) * C_num);
                        C_ad_array = C_ad_array_temp;
                        C_ad_array[C_num - 1] = id;
                        break;

                    case 8:
                        D_num++;
                        int *D_ad_array_temp = realloc(D_ad_array, sizeof(int) * D_num);
                        D_ad_array = D_ad_array_temp;
                        D_ad_array[D_num - 1] = id;
                        break;
    */
                }

            // calculate g(r) for each frame
            if(id == N - 1){
                frame++;

                // g33
                if(A_num != 0){
                    frame_AA++;
                    for(m = 0; m < A_num; m++){
                        for(n = m + 1; n < A_num; n++){
                            // compute rmn
                            r = arclength(ra[A_ad_array[m]], ra[A_ad_array[n]]);
                            // compute the frequency by checking each bin
                            for(ri = 0; ri < length_of_gr; ri++){
                                if(r <= (ri + 1) * al_interval){
                                    freq_buff33[ri]++;
                                    break;
                                }
                                else 
                                    continue;
                            }
                        }
                    }
                    for(ri = 0; ri < length_of_gr; ri++){
                        if(A_num > 1){
                            g33[ri] += (al_normalization(freq_buff33[ri], ri + 1) / (A_num * (A_num - 1) / 2));
                            freq33[ri] += (freq_buff33[ri] / (A_num * (A_num - 1) / 2));
                        }
                        else{
                            g33[ri] += al_normalization(freq_buff33[ri], ri + 1);
                            freq33[ri] += freq_buff33[ri];
                        }
                    }
                }

                // g44
                if(B_num != 0){
                    frame_BB++;
                    for(m = 0; m < B_num; m++){
                        for(n = m + 1; n < B_num; n++){
                            // compute rmn
                            r = arclength(ra[B_ad_array[m]], ra[B_ad_array[n]]);	
                            // compute the frequency by checking each bin
                            for(ri = 1; ri <= length_of_gr; ri++){
                                if(r < (ri + 1) * al_interval){
                                    freq_buff44[ri]++;
                                    break;
                                }
                                else
                                    continue;
                            }
                        }
                    }
                    for(ri = 0; ri < length_of_gr; ri++){
                        if(B_num > 1){
                            g44[ri] += (al_normalization(freq_buff44[ri], ri + 1) / (B_num * (B_num - 1) / 2));
                            freq44[ri] += (freq_buff44[ri] / (B_num * (B_num - 1) / 2));
                        }
                        else{
                            g44[ri] += al_normalization(freq_buff44[ri], ri + 1);
                            freq44[ri] += freq_buff44[ri];
                        }
                    }
                }
                    
                // g34
                if(A_num != 0 && B_num !=0){
                    frame_AB++;
                    for(m = 0; m < A_num; m++){
                        for(n = 0; n < B_num; n++){
                            // compute rmn
                            r = arclength(ra[A_ad_array[m]], ra[B_ad_array[n]]);	
                            // compute the frequency by checking each bin
                            for(ri = 1; ri <= length_of_gr; ri++){
                                if(r < (ri + 1) * al_interval){
                                    freq_buff34[ri]++;
                                    break;
                                }
                                else
                                    continue;
                            }
                        }
                    }	
                    for(ri = 0; ri < length_of_gr; ri++){
                        g34[ri] += al_normalization(freq_buff34[ri], ri + 1) / (A_num * B_num);
                        freq34[ri] += freq_buff34[ri] / (A_num * B_num);
                    }
                }

                    
            // 
            A_num = 0;
            B_num = 0;
            memset(freq_buff33, 0, sizeof(float) * length_of_gr);
            memset(freq_buff44, 0, sizeof(float) * length_of_gr);
            memset(freq_buff34, 0, sizeof(float) * length_of_gr);
            }
        }
    }

    
    // print the header of g(r) files
    fprintf(g33_out, "r g(r)\n");
    fprintf(g44_out, "r g(r)\n");
    fprintf(g34_out, "r g(r)\n");

    // normalize g(r)
    for(ri = 0; ri < length_of_gr; ri++){
        g33[ri] /= frame_AA;
        g44[ri] /= frame_BB;
        g34[ri] /= frame_AB;
        freq33[ri] /= frame_AA;
        freq44[ri] /= frame_BB;
        freq34[ri] /= frame_AB;

        fprintf(g33_out, "%f %f\n", (ri + 1) * al_interval, g33[ri]);
        fprintf(g44_out, "%f %f\n", (ri + 1) * al_interval, g44[ri]);
        fprintf(g34_out, "%f %f\n", (ri + 1) * al_interval, g34[ri]);
        fprintf(freq33_out, "%f %f\n", (ri + 1) * al_interval, freq33[ri]);
        fprintf(freq44_out, "%f %f\n", (ri + 1) * al_interval, freq44[ri]);
        fprintf(freq34_out, "%f %f\n", (ri + 1) * al_interval, freq34[ri]);
    }


    free(A_ad_array);
    free(B_ad_array);

    fclose(traj);
    fclose(g33_out);
    fclose(g44_out);
    fclose(g34_out);
    fclose(freq33_out);
    fclose(freq44_out);
    fclose(freq34_out);

/*
    fclose(g77_out);
    fclose(g88_out);
    fclose(g78_out);
    fclose(freq77_out);
    fclose(freq88_out);
    fclose(freq78_out);
*/
    return 0;
} 
