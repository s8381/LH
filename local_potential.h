// External potential \phi(z) = Az acting on the ideal gases
float local_potential(){
    extern float Lzmax;
    extern float A;
    float force = -A;

    return force;
}
