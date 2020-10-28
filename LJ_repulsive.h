// Lennard-Jones repulsive potential
float LJ_re(float sep, float epsilon_pp, int label){
    float LJ_force;
    float Rs = 2;
    float sig_A = 0.6346;
    float sig_B = 0.6441;
    float sig;
    // determine the size of the particle
    if(label == 1)
        sig = sig_A;
    else
        sig = sig_B;

    if(sep > 1 + Rs) // cut-off = 1 nm
        LJ_force = 0;
    else {
        float s = (sig / (sep - Rs)) * (sig / (sep - Rs)) * (sig / (sep - Rs));
        float S = s * s;
        LJ_force = 48 * epsilon_pp * S * S / ((sep - Rs) * sep);
    }
    return(LJ_force);
}
