#include "matrice_GS.hpp"

//---------------------------------------------------
//             Creations, destructions
//---------------------------------------------------

//------------ createur ------------
mat::mat()
{
  row = 0;
  col = 0;
  GScoord = NULL;
}

//------------ destructeur ------------
mat::~mat()
{
  if (GScoord!=NULL) 
    delete[] GScoord;
  GScoord = NULL;
  row = 0;
  col = 0;
}


// modif Fred: operateur de recopie
mat::mat(const mat & M)
{
  row = M.row;
  col = M.col;
  if (M.GScoord!=NULL) {
	GScoord = new double[col*row];
	for (int i=0;i<(row*col);i++) {
	  GScoord[i] = M.GScoord[i];
	}
  } else GScoord=NULL;
}


//------------ initialisation ------------
void mat::set_size(int Row, int Col)
  // met le mat courant à la longueur souhaitée
{
  row = Row;
  col = Col;
  if (GScoord != NULL) 
    delete[] GScoord;
  GScoord = new double[col*row];
}

mat::mat(int Row, int Col)
  //construit un mat de taille donnée
{
  row = Row;
  col = Col;
  GScoord = new double[Row*Col];
}

//------------ creation meme matrice ------------
mat & mat::operator= (const mat & a)
// affectation surchargée
{
  if (this!=&a) {
    row = a.row;
    col = a.col;
    if (GScoord != NULL) 
      delete[] GScoord;
    if (a.GScoord != NULL) 
      {
	GScoord = new double[row*col];
	for (int i=0;i<row*col;i++) 
	  GScoord[i] = a.GScoord[i];
      } else {
	GScoord = NULL;
      }
  }
  return *this;
}

//------------ mise a zero ------------
void mat::zeros()
  // met le mat courant à 0
{
  if (GScoord!=NULL) {
    for (int i=0; i<col*row; i++) 
      GScoord[i]=0.0;
  }
}

//------------ acceder à un élément ------------
double & mat::operator() (int a,int b)
  // renvoie la composante d'un mat
{
  if (GScoord != NULL) {
    return GScoord[a*col+b];
  } else {
    cerr<<"on demande la composante d'un mat vide"<<endl;
    exit(1);
  }
}

//---------------------------------------------------
//         Opérations élémentaires surchargés
//---------------------------------------------------

//------------ addition de deux matrices ------------
mat mat::operator+ (const mat & a)
  // addition surchargée
{
  mat ret;
  if ( (row != a.row) && (col != a.col) )
    {
      cerr<<"addition de mats de tailles differentes impossible"<<endl;
      exit(1);
    } else {
      if ((GScoord != NULL) && (a.GScoord != NULL)) {
	ret.row = a.row;
	ret.col = a.col;
	ret.GScoord = new double[col*row];
	for (int i=0; i<col*row; i++) 
	  ret.GScoord[i] = GScoord[i] + a.GScoord[i];
      }
    }
  return ret;
}

//------------ soustraction ------------
mat mat::operator- (const mat & a)
  // addition surchargée
{
  mat ret;
  if ( (row != a.row) && (col != a.col) )
    {
      return ret;
    } else {
      if ((GScoord != NULL) && (a.GScoord != NULL)) {
	ret.row = a.row;
	ret.col = a.col;
	ret.GScoord = new double[col*row];
	for (int i=0; i<col*row; i++) 
	  ret.GScoord[i] = GScoord[i] - a.GScoord[i];
      }
    }
  return ret;
}

//------------ multiplication par une constante ------------
mat mat::operator* (const double c)
  // multiplication par un scalaire
{
  mat ret;
  if (GScoord != NULL) {
    ret.row = row;
    ret.col = col;
    ret.GScoord = new double[col*row];
    for (int i=0; i<col*row; i++) 
      ret.GScoord[i] = GScoord[i]*c;
  }
  return ret;
}

//------------ multiplication à droite par a ------------
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

// Modif Fred: matrice par vecteur
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

// Modif Fred: la trace
double mat::trace()
{
  double tr = 0.;
  
  if (GScoord != NULL) {
	if (row != col) {
	  cerr<<"on demande la trace d'un mat non carre"<<endl;
	  exit(1);
	} else {
	  for (int i=0; i<col; i++) {
		tr += GScoord[i*col+i];
	  }
	  return tr;
	}
  } else {
    cerr<<"on demande la trace d'un mat vide"<<endl;
    exit(1);
  }
}


