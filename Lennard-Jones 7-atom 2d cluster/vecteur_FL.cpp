#include "vecteur_FL.hpp"

/////////////////////////////////////
// METHODS of the vector class


vec::vec()
  // default constructor
{
  //  cerr<<"simple creation"<<endl;
  FLdim = 0;
  FLcoord = NULL;
}


vec::vec(int d)
  //create a vector of given size 'd'
{
  //  cerr<<"deatiled creation"<<d<<endl;
  FLdim = d;
  FLcoord = new double[d];
}


vec::vec(const vec & v)
  // copy 'v' into current
{
  //  cerr<<"copy"<<endl;
  FLdim = v.FLdim;
  if (v.FLcoord!=NULL)
  {
	FLcoord = new double[FLdim];
	for (int i=0;i<FLdim;i++) FLcoord[i]=(v.FLcoord[i]);
  }
  else FLcoord=NULL;
}


vec::~vec()
  // destructor
{
  //  cerr<<"destructor"<<endl;
  if (FLcoord!=NULL) delete[] FLcoord;
  FLcoord = NULL;
  FLdim =0;
}


vec & vec::operator= (const vec & a)
// overload assignment for =
{
  if (this!=&a) {
	//	cerr<<"do something ="<<endl;
	FLdim = a.FLdim;
	if (FLcoord!=NULL) delete[] FLcoord;
	if (a.FLcoord!=NULL)
    {
	  FLcoord = new double[FLdim];
	  for (int i=0;i<FLdim;i++) FLcoord[i]=a.FLcoord[i];
	}
    else
    {
	  FLcoord = NULL;
	}
  }
  else
  {
	//	cerr<<"do not do this ="<<endl;
  }
  return *this;
}


vec vec::operator+ (const vec & a)
  // overload +
{
  //  cerr<<"enter in +"<<endl;

  vec b;

  if (FLdim!=a.FLdim)
  {
	cerr<<"addition of vectors of different sizes is impossible"<<endl;
	exit(1);
  }
  else
  {
	if ((FLcoord != NULL) && (a.FLcoord!=NULL))
    {
	  b.FLdim = a.FLdim;
	  b.FLcoord = new double[FLdim];
	  for (int i=0;i<FLdim;i++) b.FLcoord[i] = FLcoord[i] + a.FLcoord[i];
	}
  }
  return b;
}

vec vec::operator/ (const double a)
  // division by a scalar
{
  vec b;

  if (FLcoord != NULL)
  {
    b.FLdim = FLdim;
    b.FLcoord = new double[FLdim];
    for (int i=0;i<FLdim;i++) b.FLcoord[i] = FLcoord[i]/a;
  }

  return b;
}

double & vec::operator() (int a)
  // return a component of the vector
{
  if (FLdim-1<a)
  {
	cerr<<"this component does not exist - vector issue"<<endl;
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
	  cerr<<"this vector is null"<<endl;
	  exit(1);
	}
  }
}

double & vec::operator[] (int a)
  // same as ()
{
  if (FLdim-1<a)
  {
	cerr<<"this component does not exist "<<endl;
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
	  cerr<<"this vector is null"<<endl;
	  exit(1);
	}
  }
}


void vec::set_size(int d)
  // set current vector to the desired length
{
  FLdim = d;
  if (FLcoord!=NULL) delete[] FLcoord;
  FLcoord = new double[FLdim];
}

void vec::zeros()
  // make a vector 0
{
  if (FLcoord!=NULL)
  {
	for (int i =0;i<FLdim;i++) FLcoord[i]=0.0;
  }
}


// vec vec::mult(double a)
//   // multiplication of the current vector by a scalar and return the result
// {
//   vec b(FLdim);
//   if (FLcoord!=NULL)
//   {
// 	    for (int i =0;i<FLdim;i++) b.FLcoord[i]=a*(FLcoord[i]);
//   }
//   return b;
// }

double vec::scal(const vec & v)
  // scalar product of two vectors
{
  double p = 0;
  if (FLdim!=v.FLdim)
  {
	cerr<<"Scalar product requires the vectors to be of same dimensions"<<endl;
	exit(1);
  }
  else
  {
	if ((FLcoord!=NULL) && (v.FLcoord!=NULL))
    {
	  for (int i=0;i<FLdim;i++)
      {
		p += (FLcoord[i])*(v.FLcoord[i]);
	  }
	}
  }
  return p;
}

double vec::norm2()
  // L2 norm of a vector
{
  double p = 0.;
  for (int i=0;i<FLdim;i++)
  {
	p += (FLcoord[i])*(FLcoord[i]);
  }
  return p;
}

// Additional functions - scalar multiplication by a vector

vec operator*(double t, const vec &v)
{
    int i;
    vec r(v.FLdim);

    for (i=0; i<v.FLdim; i++)
	r.FLcoord[i] = t * v.FLcoord[i];

    return r;
}
