float LJ_shift(float sep){
    extern float Rs;
    extern float epsilon;
    extern float sig;
    float LJ_force;


    if(sep > 2.5 * sig + Rs)
        LJ_force = 0;
    else{
        float s = (sig / (sep - Rs)) * (sig / (sep - Rs)) * (sig / (sep - Rs));
	float S =  s * s; 
        LJ_force = 24 * epsilon * (2 * S * S - S) / ((sep - Rs) * sep);
    }

    return(LJ_force);
}
