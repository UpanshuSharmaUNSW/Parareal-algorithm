#include "gaussien.hpp"
#include <math.h>
#include <stdlib.h>

gaussien::gaussien()
{
    random_flag = 1;
    //mtrand1.seed(1);
    mtrand1.seed(1);
}

// le generateur uniforme entre 0 et 1
double gaussien :: random_number()
{
    //double res = ((double)rand()/RAND_MAX);
    //  double res = drand48();
    double res = mtrand1.randDblExc();
    return res;
}

// le generateur de va gaussienne de variance "SIGMA^2": cette methode suit la methode polaire,
double gaussien :: rand_number(double m, double sigma){
    double temp1,temp2;

    if (random_flag){
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

// Methode de gaussienne tronquee = tirage sur [-dm,dm]
// cet algorithme est fonde sur une methode de rejet...
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
