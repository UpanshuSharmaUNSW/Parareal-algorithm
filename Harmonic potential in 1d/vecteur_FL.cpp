#include "vecteur_FL.hpp"

/////////////////////////////////////
// METHODES de la classe vecteur


vec::vec()
// constructeur par defaut
{
    //  cerr<<"creation simple"<<endl;
    FLdim = 0;
    FLcoord = NULL;
}


vec::vec(int d)
//construit un vec de taille donn�e
{
    //  cerr<<"creation sophistiqu� "<<d<<endl;
    FLdim = d;
    FLcoord = new double[d];
}


vec::vec(const vec & v)
// operateur de recopie de v dans le courant
{
    //  cerr<<"recopie"<<endl;
    FLdim = v.FLdim;
    if (v.FLcoord!=NULL) {
        FLcoord = new double[FLdim];
        for (int i=0;i<FLdim;i++) FLcoord[i]=(v.FLcoord[i]);
    } else FLcoord=NULL;
}


vec::~vec()
// destruction
{
    //  cerr<<"destruction"<<endl;
    if (FLcoord!=NULL) delete[] FLcoord;
    FLcoord = NULL;
    FLdim =0;

}


vec & vec::operator= (const vec & a)
// affectation surcharg�e
{
    if (this!=&a) {
        //	cerr<<"faire qqch dans ="<<endl;
        FLdim = a.FLdim;
        if (FLcoord!=NULL) delete[] FLcoord;
        if (a.FLcoord!=NULL) {
            FLcoord = new double[FLdim];
            for (int i=0;i<FLdim;i++) FLcoord[i]=a.FLcoord[i];
        } else {
            FLcoord = NULL;
        }
    } else {
        //	cerr<<"ne rien faire dans ="<<endl;
    }
    return *this;
}


vec vec::operator+ (const vec & a)
// addition surcharg�e
{
    //  cerr<<"entree dans +"<<endl;

    vec b;

    if (FLdim!=a.FLdim) {
        cerr<<"addition de vecs de tailles differentes impossible"<<endl;
        exit(1);
    } else {
        if ((FLcoord != NULL) && (a.FLcoord!=NULL)) {
            b.FLdim = a.FLdim;
            b.FLcoord = new double[FLdim];
            for (int i=0;i<FLdim;i++) b.FLcoord[i] = FLcoord[i] + a.FLcoord[i];
        }
    }
    return b;
}

vec vec::operator/ (const double a)
// division par un scalaire
{
    vec b;

    if (FLcoord != NULL) {
        b.FLdim = FLdim;
        b.FLcoord = new double[FLdim];
        for (int i=0;i<FLdim;i++) b.FLcoord[i] = FLcoord[i]/a;
    }

    return b;
}

// we use this routine in the code
double & vec::operator() (int a)
// renvoie la composante d'un vec
{
    if (FLdim-1<a)
    {
        cerr<<"Requested component of vector does not exist"<<endl;
        cerr<<"Requested component #=" << a << ", and size of vector = " << FLdim-1 <<endl;
        exit(1);
    }
    else
    {
        if (FLcoord!=NULL)
        {
            return FLcoord[a];
        }
        else
        {
            cerr<<"Component requested from empty vector"<<endl;
            exit(1);
        }
    }
}

double & vec::operator[] (int a)
// idem que ()
{

    if (FLdim-1<a)
    {
        cerr<<"Requested component of vector does not exist"<<endl;
        cerr<<"Requested component #=" << a << ", and size of vector = " << FLdim-1 <<endl;
        exit(1);
    }
    else
    {
        if (FLcoord!=NULL)
        {
            return FLcoord[a];
        }
        else
        {
            cerr<<"Component requested from empty vector"<<endl;
            exit(1);
        }
    }
}


void vec::set_size(int d)
// met le vec courant � la longueur souhait�e
{
    FLdim = d;
    if (FLcoord!=NULL) delete[] FLcoord;
    FLcoord = new double[FLdim];
}

void vec::zeros()
// met le vec courant � 0
{
    if (FLcoord!=NULL) {
        for (int i =0;i<FLdim;i++) FLcoord[i]=0.0;
    }
}


// vec vec::mult(double a)
//   // multiplie le vec courant par un scalaire et renvoie le resultat
// {
//   vec b(FLdim);
//   if (FLcoord!=NULL) {
// 	for (int i =0;i<FLdim;i++) b.FLcoord[i]=a*(FLcoord[i]);
//   }
//   return b;
// }

double vec::scal(const vec & v)
// defintion du produit scalaire
{
    double p = 0;
    if (FLdim!=v.FLdim) {
        cerr<<"vecs de FLdim diff dans scal"<<endl;
        exit(1);
    } else {
        if ((FLcoord!=NULL) && (v.FLcoord!=NULL)) {
            for (int i=0;i<FLdim;i++) {
                p += (FLcoord[i])*(v.FLcoord[i]);
            }
        }
    }
    return p;
}

double vec::norm2()
// defintion du produit scalaire
{
    double p = 0.;
    for (int i=0;i<FLdim;i++) {
        p += (FLcoord[i])*(FLcoord[i]);
    }

    return p;
}

// Fonctions supplementaires

vec operator*(double t, const vec &v)
{
    int i;
    vec r(v.FLdim);

    for (i=0; i<v.FLdim; i++)
    r.FLcoord[i] = t * v.FLcoord[i];

    return r;
}
