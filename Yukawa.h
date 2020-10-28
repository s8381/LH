// Yukawa potential
float yukawa(float sep, float miu, int label){
    extern float Rs;
    extern float sig_A;
    extern float sig_B;
    extern float kappa;
    float sig;
    float force;

    if(label == 1)
        sig = sig_A;
    else
        sig = sig_B;

    if(sep > 6 / kappa)
        force = 0;
    else
        force = miu * (Rs + sig / 2) * exp(-kappa * (sep - (Rs + sig / 2))) * (kappa + 1 / sep) / (sep * sep); 

    return force;
}
