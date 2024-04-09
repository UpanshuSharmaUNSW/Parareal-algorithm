#include "gaussien.hpp"
#include <math.h>
#include <stdlib.h>

gaussien::gaussien()
{
  random_flag = 1;
  //mtrand1.seed(1);
  mtrand1.seed(1);
  // mtrand1.seed(1234); // Use different seed for gamma=10
}

 // uniform distribution between 0 and 1
double gaussien :: random_number()
{

  //double res = ((double)rand()/RAND_MAX);
  //  double res = drand48();
  double res = mtrand1.randDblExc();
  return res;

}

 // Gaussian variance generator of "SIGMA^2": this method uses the polar method
double gaussien :: rand_number(double m, double sigma)
{
double temp1,temp2;

 if (random_flag)
 {
   temp1=random_number();
   temp2=random_number();
   if (temp1==0) temp1=1e-9;
   r1=m+sigma*(sqrt(-2.*log(temp1))*cos(2.*M_PI*temp2));
   r2=m+sigma*(sqrt(-2.*log(temp1))*sin(2.*M_PI*temp2));
   random_flag=0;
   return r1;
 }
 else{
   random_flag=1;
   return r2;
 }
}

// Truncated Gaussian method based on rejection method = drawing from [-dm,dm]
double gaussien::rand_number_tronque(double m, double sigma, double dm)
{
  bool echec = true;
  double res;
  while (echec)
  {
    res = rand_number(m,sigma);
    if ( (res < dm) && (res > -dm) )
	{
	  echec = false;
	}
  }

  return res;

}
