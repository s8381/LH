// Dot product function
float dot_product(float  *a,float *b) {
    float result = 0;
    int p;
    for(p = 0; p < 3; p++)
        result += *(a + p) * *(b + p);
    return result;
}