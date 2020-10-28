// Debye-Hueckel potential
float DH(float sep, int label_A, int label_B){
    // cutoff = 3
    extern float qa;
    extern float qb;
    extern float sig_A;
    extern float sig_B;
    extern float kappa;
    extern float lb;
    extern float rcut;
    float force;
    float phi0;
    float sig;

    if(label_A == 1 && label_B == 1){
        phi0 = qa * qa * lb;
        sig = sig_A;
    }
    else if(label_A == 2 && label_B == 2){
        phi0 = qb * qb * lb;
        sig = sig_B;
    }
    else{
        phi0 = qa * qb * lb;
        sig = (sig_A + sig_B) / 2;
    }

    if(sep > rcut)
        force = 0;
    else
        force = phi0 * exp(-kappa * (sep - sig)) * (1 / sep + kappa) / (sep * sep * (1 + kappa * sig / 2) * (1 + kappa * sig / 2));

    return(force);
}
