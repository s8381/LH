// wall potential
float wall(float z){
    extern float epsilon_pp;
    extern float sig0;
    extern float Lzmax;
    float LJ_force, sep;

    if(z <= 0)
        sep = z + (Lzmax + sig0 / 2);
    else
        sep = z - (Lzmax + sig0 / 2);

    if(sep > 2.5 * sig0 || sep < -2.5 * sig0) // cut-off = 2.5 * sig
        LJ_force = 0;
    else {
	float S = (sig0 / sep) * (sig0 / sep) * (sig0 / sep) * (sig0 / sep) * (sig0 / sep) * (sig0 / sep);
        LJ_force = 24 * epsilon_pp * (2 * S * S) / sep;
    }

    return(LJ_force);
}
