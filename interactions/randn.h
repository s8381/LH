// Normal distribution random number
float randn (float mu, float sigma){
    float U1, U2, W, mult;
    static float X1, X2;
    static float call = 0;
 
    if (call == 1){
        call = !call;
        return (mu + sigma * (float) X2);
    }
 
  do
    {
      U1 = -1 + ((float) rand () / RAND_MAX) * 2;
      U2 = -1 + ((float) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;

  return (mu + sigma * (float) X1);
}
