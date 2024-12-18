/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{
  if (*kv){++*lab;}
  {
    AB[0,0]=0;
    AB[*kv,0]=0;
    AB[*kv+1,0]=2;
    AB[*kv+2,0]=-1;
    for (int i=1; i<*la-1; ++i)
    {
      AB[0,i]=0;
      AB[*kv,i]=-1;
      AB[*kv+1,i]=2;
      AB[*kv+2,i]=-1;
    }
    AB[0,*la]=0;
    AB[*kv,*la]=-1;
    AB[*kv+1,*la]=2;
    AB[*kv+2,*la]=0;
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv)
{
  if (*kv){++*lab;}
  {
    AB[0,0]=0;
    AB[*kv,0]=0;
    AB[*kv+1,0]=1;
    AB[*kv+2,0]=0;
    for (int i=1; i<*la-1; ++i)
    {
      AB[0,i]=0;
      AB[*kv,i]=0;
      AB[*kv+1,i]=1;
      AB[*kv+2,i]=0;
    }
    AB[0,*la]=0;
    AB[*kv,*la]=0;
    AB[*kv+1,*la]=1;
    AB[*kv+2,*la]=0;
  }
}

//Vecteur initial (T0 en première valeur et T1 à la fin, avec que des zéros au milieu)
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
  RHS[0]=*BC0;
  for (int i=1; i<*la-1; ++i)
  {
    RHS[i]=0;
  }
  RHS[0]=*BC1;
}  

//Remplir un vecteur avec les valeurs qui seraient données par la solution analytique (T(x)=T0 + x(T1 - T0))
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1)
{
  int n = *la;         
  double T0 = *BC0;     
  double T1 = *BC1;     
  
  for (int i = 0; i < n; i++) 
  {
  EX_SOL[i] = T0 + X[i] * (T1 - T0);
  }
}  

void set_grid_points_1D(double* x, int* la)
{
 int n = *la;  
 double dx = 1.0 / (n - 1); 

  for (int i = 0; i < n; i++) 
  {
    x[i] = i * dx;
  }
}

double relative_forward_error(double* x, double* y, int* la)
{
  int n = *la;  
  double norm_diff = 0.0; 
  double norm_x = 0.0;    
  
  for (int i = 0; i < n; i++) 
  {
    double diff = x[i] - y[i];
    norm_diff += diff * diff; 
    norm_x += x[i] * x[i];    
  }

  norm_diff = sqrt(norm_diff); 
  norm_x = sqrt(norm_x);       
  
  if (norm_x == 0.0) 
  {
    return -1.0; 
  }
  return norm_diff / norm_x;
}

int indexABCol(int i, int j, int *lab)
{
  int kl = (*lab - 1) / 2; 
  int localRow = kl + i - j; 

  if (localRow < 0 || localRow >= *lab) 
  {
    return -1;
  }

  return localRow + j * (*lab);
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  int i, j, k;
    double factor;
    *info = 0;    //Veut dire qu'il n'y a pas d'erreur.

    //Factorisation LU
    for (i = 0; i < *n; i++) 
    {
      //MAJ de la diagonale principale de A (AB[0])
      if (i > 0) 
      {
        AB[0 + i * (*lab)] -= AB[1 + (i - 1) * (*lab)] * AB[2 + (i - 1) * (*lab)] / AB[0 + (i - 1) * (*lab)];
      }

      //Calcul du pivot et mise à jour des éléments dans les bandes
      if (i < *n - 1) 
      {
        //Calcul de la factorisation de la diagonale inf (AB[1])
        AB[1 + i * (*lab)] = AB[1 + i * (*lab)] / AB[0 + i * (*lab)];
        //MAJ de la diagonale sup (AB[2])
        AB[2 + (i + 1) * (*lab)] -= AB[1 + i * (*lab)] * AB[2 + i * (*lab)];
      }
      ipiv[i] = i + 1; //Pas de permutation ici, indice inchangé
    }

    //Vérification de l'inversibilité
    for (i = 0; i < *n; i++) 
    {
      if (AB[0 + i * (*lab)] == 0) 
      {
        *info = i + 1; //Indique l'indice où l'élément pivot est nul
        return 0;
      }
    }
  *info = 0;
  return *info;
}