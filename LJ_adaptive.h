// Lennard-Jones potential
float LJ(float sep, float za, float zb){
    extern float epsilon_pp;
    extern float f_scaling;
    extern float Lzmax;
    extern float sig0;
    float LJ_force;
    float sig, sig_A, sig_B;

    sig_A = sig0 + f_scaling * (za + Lzmax) / (2 * Lzmax) * sig0;
    sig_B = sig0 + f_scaling * (zb + Lzmax) / (2 * Lzmax) * sig0;
    sig = (sig_A + sig_B) / 2;

    if(sep > 2.5 * sig) // cut-off = 2.5 * sig
        LJ_force = 0;
    else {
        float s = (sig / sep) * (sig / sep) * (sig / sep);
	float S = s * s;
        LJ_force = 24 * epsilon_pp * (2 * S * S - S) / (sep * sep);
    }
    return(LJ_force);
}
