// 9-3 Lennard-Jones potential
float LJ93(float sep, int label_a, int label_b){
    extern float epsilon_pp;
    extern float sig_A;
    extern float sig_B;
    float LJ_force;
    float sig;

    if(label_a == 1 && label_b ==1)
        sig = sig_A;
    else if(label_a == 2 && label_b == 2)
        sig = sig_B;
    else
        sig = (sig_A + sig_B) / 2;

    if(sep > 2.5 * sig) // cut-off = 2.5 * sig
        LJ_force = 0;
    else{
        float s = (sig / sep) * (sig / sep) * (sig / sep);
        float S = s * s * s;
        LJ_force = 3 * sqrt(3) * epsilon_pp * (9 * S - 3 * s) / (2 * sep * sep);
    }
    return(LJ_force);
}
