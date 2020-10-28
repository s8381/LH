// Lennard-Jones potential
double LJ6(double sep, double epsilon_pp){
    double LJ_force;
    double sig = 0.35;
    if(sep > 2.5 * sig) // cut-off = 2.5 sigma
        LJ_force = 0;
    else {
        double s = (sig / sep) * (sig / sep) * (sig / sep);
	double S = s * s;
        LJ_force = 24 * epsilon_pp * S  / (sep * sep);
    }
    return(LJ_force);
}