// tengential harmonic potential
float *t_harmonic(float *u_vec, float *ra){
    extern float kh;
    int i;
    float t_dis[3];
    float *force = malloc(sizeof(float) * 3);
    float dot = 0;

    // ra fot n
    for(i = 0; i < 3; i++)    
        dot += ra[i] * u_vec[i];

    // rt = ra - rn = ra - (ra dot n)n
    for(i = 0; i < 3; i++){
        t_dis[i] = ra[i] - dot * u_vec[i];
    }

    for(i = 0; i < 3; i++){
        force[i] = - kh * t_dis[i];
    }

    return force;
}
