#include "matrice_GS.hpp"

//---------------------------------------------------
//             Constructors, Destructors
//---------------------------------------------------

//------------ Constructor ------------
mat::mat()
{
  row = 0;
  col = 0;
  GScoord = NULL;
}

//------------ Destructor ------------
mat::~mat()
{
  if (GScoord!=NULL)
    delete[] GScoord;
  GScoord = NULL;
  row = 0;
  col = 0;
}


// Copy operator
mat::mat(const mat & M)
{
  row = M.row;
  col = M.col;
  if (M.GScoord!=NULL)
  {
	GScoord = new double[col*row];
	for (int i=0;i<(row*col);i++)
    {
	  GScoord[i] = M.GScoord[i];
	}
  } else GScoord=NULL;
}


//------------ initialisation ------------
void mat::set_size(int Row, int Col)
  // set current matrix to desired dimensions
{
  row = Row;
  col = Col;
  if (GScoord != NULL)
    delete[] GScoord;
  GScoord = new double[col*row];
}

mat::mat(int Row, int Col)
  // construct a matrix of given size
{
  row = Row;
  col = Col;
  GScoord = new double[Row*Col];
}

//------------ creation of matrix ------------
mat & mat::operator= (const mat & a)
// overloaded assignment
{
  if (this!=&a)
  {
    row = a.row;
    col = a.col;
    if (GScoord != NULL)
      delete[] GScoord;
    if (a.GScoord != NULL)
    {
        GScoord = new double[row*col];
        for (int i=0;i<row*col;i++)
        GScoord[i] = a.GScoord[i];
    }
    else
    {
        GScoord = NULL;
    }
  }
  return *this;
}

//------------ mise a zero ------------
void mat::zeros()
  // set current matrix to 0
{
  if (GScoord!=NULL)
  {
    for (int i=0; i<col*row; i++)
      GScoord[i]=0.0;
  }
}

//------------ access an element ------------
double & mat::operator() (int a,int b)
  // returns the specified component of matrix
{
  if (row-1<a)
    {
      cerr<<"this component " <<a<<" does not exist, larger than row "<<row << "(matrix problem)" << endl;
      exit(1);
    }
  else
    {
      if (col-1<b)
	{
	  cerr<<"this component "<<b<<" does not exist, larger than col "<<col<< "(matrix problem)" << endl;
	  exit(1);
	}
      else
	{
	  if (GScoord != NULL)
	    {
	      return GScoord[a*col+b];
	    }
	  else
	    {
	      cerr<<"The matrix is empty - no element returned"<<endl;
	      exit(1);
	    }
	}
    }
}

//---------------------------------------------------
//         Overloaded elementary operations
//---------------------------------------------------

//------------ addition of two matrices ------------
mat mat::operator+ (const mat & a)
  // addition overloaded
{
  mat ret;
  if ( (row != a.row) && (col != a.col) )
    {
      cerr<<"matrices are of different dimensions and cannot be added"<<endl;
      exit(1);
    }
  else {
      if ((GScoord != NULL) && (a.GScoord != NULL))
      {
          ret.row = a.row;
          ret.col = a.col;
          ret.GScoord = new double[col*row];
          for (int i=0; i<col*row; i++)
              ret.GScoord[i] = GScoord[i] + a.GScoord[i];
      }
    }
  return ret;
}

//------------ substraction ------------
mat mat::operator- (const mat & a)
  // addition surcharg�e
{
  mat ret;
  if ( (row != a.row) && (col != a.col) )
  {
      return ret;
  }
  else
  {
      if ((GScoord != NULL) && (a.GScoord != NULL))
      {
          ret.row = a.row;
          ret.col = a.col;
          ret.GScoord = new double[col*row];
          for (int i=0; i<col*row; i++)
              ret.GScoord[i] = GScoord[i] - a.GScoord[i];
      }
   }
  return ret;
}

//------------ multiplication by a constant ------------
mat mat::operator* (const double c)
  // multiplication by a scalaire
{
  mat ret;
  if (GScoord != NULL)
  {
    ret.row = row;
    ret.col = col;
    ret.GScoord = new double[col*row];
    for (int i=0; i<col*row; i++)
      ret.GScoord[i] = GScoord[i]*c;
  }
  return ret;
}

//------------ multiplication from the right by a matrix ------------
mat mat::operator*(const mat & a)
{
  mat ret;
  if(a.row != col)
    return ret;
  ret.set_size(row,a.col);
  ret.zeros();
  for(int i = 0; i < row; i++)
    for(int j = 0; j < a.col; j++)
      for(int k = 0; k < col; k++)
          ret(i,j) = ret(i,j) + GScoord[i*col+k]*a.GScoord[k*a.col+j];
  return ret;
}

// Multiplication with a vector from right
vec mat::mult(const vec & a)
{
  vec ret;
  if(a.FLdim != col)
    return ret;
  ret.set_size(row);
  ret.zeros();
  for(int i = 0; i < row; i++)
	for(int k = 0; k < col; k++)
	  ret(i) = ret(i) + GScoord[i*col+k]*a.FLcoord[k];
  return ret;
}

// Trace of matrix  - sum of diagonal element
double mat::trace()
{
  double tr = 0.;

  if (GScoord != NULL)
  {
	if (row != col)
    {
	  cerr<<"trace of non-square matrix not defined"<<endl;
	  exit(1);
	}
    else
    {
	  for (int i=0; i<col; i++)
      {
		tr += GScoord[i*col+i];
	  }
	  return tr;
	}
  }
  else
  {
    cerr<<"trace of empty matrix not defined"<<endl;
    exit(1);
  }
}
