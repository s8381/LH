// Lennard-Jones potential
float LJ_local(float sep, float sig_A, float sig_B){
    extern float epsilon_pp;
    float LJ_force;
    float sig;

    // calculate the locality-dependent radii 
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
