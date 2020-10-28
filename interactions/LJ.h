// Lennard-Jones potential
float LJ(float sep){
    extern float epsilon_pp;
    extern float sig;
    float LJ_force;

    if(sep > 2.5 * sig) // cut-off = 2.5 * sig
        LJ_force = 0;
    else {
        float s = (sig / sep) * (sig / sep) * (sig / sep);
	float S = s * s;
        LJ_force = 24 * epsilon_pp * (2 * S * S - S) / (sep * sep);
    }
    return(LJ_force);
}
