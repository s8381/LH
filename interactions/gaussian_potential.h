// Gaussian potential
float gaussian(float sep, int label_a, int label_b){
    extern float sig_A;
    extern float sig_B;
    float gaussian_force;
    float sig;

    if(label_a == 1 && label_b ==1)
        sig = sig_A;
    else if(label_a == 2 && label_b == 2)
        sig = sig_B;
    else
        sig = (sig_A + sig_B) / 2;

    if(sep > 2.5 * sig) // cut-off = 2.5 * sig
        gaussian_force = 0;
    else
        gaussian_force = 2 * exp(-(sep / sig) * (sep / sig)) / sig / sig;

    return(gaussian_force);
}
