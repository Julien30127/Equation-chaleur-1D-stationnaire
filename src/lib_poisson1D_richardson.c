/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la)
{
  for (int k = 0; k < *la; k++) 
  {
    eigval[k] = 4.0 * pow(sin(k * M_PI / (2.0 * (*la + 1))), 2);
  }
}

double eigmax_poisson1D(int *la)
{
  double h = 1.0/(*la + 1);
  double *eigval=malloc(*la*sizeof(double));
  eig_poisson1D(eigval, *la);
  double max=0.0;
  for int(int k=0; k<*la; k++)
  {
    max = fmax(eigval[k], max);
  }
  free (eigval); 
  return max;
}

double eigmin_poisson1D(int *la)
{
  double h = 1.0/(*la + 1);
  double *eigval;
  eig_poisson1D(eigval, *la);
  double min=0.0;
  for int(int k=0; k<*la; k++)
  {
    min = fmin(eigval[k], min);
  }  
  free (eigval);
  return min;
}

double richardson_alpha_opt(int *la)
{
  double lambda_min = eigmin_poisson1D(*la);
  double lambda_max = eigmax_poisson1D(*la);
  double alpha = 2/(lambda_min+lambda_max);
  double x_k=0.0;
  double x_kp1=0.0;


  return 0;
}

/*EXERCICE 7*/

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
  //Allocation du résidu
  double *y = malloc(*la * sizeof(double));  
  //Norme du résidu
  double res_norm;
  //Norme du vecteur RHS                          
  double rhs_norm;                          

  //Initialisation

  cblas_dcopy(*la, RHS, 1, y, 1);           //Equivalent à y = RHS
  rhs_norm = cblas_dnrm2(*la, RHS, 1);      //Norme initiale RHS
  if (rhs_norm == 0.0) {rhs_norm = 1.0;}      //Si valeur interdite
  *nbite = 0;

  //Ouverture du fichier où seront stockées les valeurs du résidu

  FILE *file = fopen("residus.txt", "w");
  if (file == NULL) 
  {
    perror("Erreur lors de l'ouverture du fichier");
    free(y);
    return;
  }

  for (int k = 0; k < *maxit; k++) 
  {
    //Calcul du résidu
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, y, 1);

    //Calcul de la norme du résidu
    res_norm = cblas_dnrm2(*la, y, 1);
    resvec[k] = res_norm / rhs_norm;      
    printf("Itération %d : Résidu normalisé = %.8f\n", k, resvec[k]);

    //Ecriture dans le fichier
    fprintf(file, "%d %.8f\n", k, resvec[k]);
    printf("Itération %d : Résidu normalisé = %.8f\n", k, resvec[k]);

    //Test de convergence
    if (resvec[k] < *tol) 
    {
      *nbite = k + 1;                   
      break;
    }

    //Mise à jour de X
    cblas_daxpy(*la, *alpha_rich, y, 1, X, 1);
    *nbite = k + 1;
  }
  //Fermeture du fichier
  fclose(file);

  free(y);
}

double compute_error_Richardson(double *X, double *X_exact, int la)
{
  double *e = malloc(la * sizeof(double));
  cblas_dcopy(la, X, 1, e, 1);               //Copie X dans e
  cblas_daxpy(la, -1.0, X_exact, 1, e, 1);  //Equivalent à e = e - X_exact
  double error = cblas_dnrm2(la, e, 1);     //Norme de l'erreur
  free(e);
  return error;
}

/*EXERCICE 8*/

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv)
{
  //Initialisation à zéro
  for (int i = 0; i < (*la) * (*kv); i++) 
  {
    MB[i] = 0.0;
  }

  //Copie de la diagonale principale
  for (int i = 0; i < *la; i++) 
  {
    MB[i * (*kv) + *kl] = AB[i * (*lab) + *kl]; //La diagonale est décalée de *kl dans AB
  }
}

void jacobi_tridiag(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
  int kv = *kl + *ku + 1; // Taille effective d'une ligne dans AB
  double *MB = malloc((*la) * kv * sizeof(double));   //Matrice diagonale M
  double *BX = malloc(*la * sizeof(double));        //Vecteur temporaire BX
  double *residual = malloc(*la * sizeof(double));  //Vecteur de résidus
  double rhs_norm, res_norm;

  //On extrait la matrice diagonale principale à l'aide de la fonction précédemment implémentée.
  extract_MB_jacobi_tridiag(AB, MB, lab, la, ku, kl, &kv);

  //Initialisation

  rhs_norm = cblas_dnrm2(*la, RHS, 1);
  if (rhs_norm == 0.0) {rhs_norm = 1.0;}

  //Même principe que pour l'exercice 7

  FILE *file = fopen("résidus_jacobi.txt", "w");

  if (file == NULL) 
  {
    perror("Erreur lors de l'ouverture du fichier");
    return;
  }

  for (int k = 0; k < *maxit; k++) 
  {
    //BX = RHS - (A - M) * X
    cblas_dcopy(*la, RHS, 1, BX, 1); // BX = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, BX, 1);

    //X = M^{-1} * BX (division élément par élément)
    for (int i = 0; i < *la; i++) 
    {
      X[i] = BX[i] / MB[i * kv + *kl];
    }

    //Calcul du résidu

    cblas_dcopy(*la, RHS, 1, residual, 1); // residual = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, residual, 1);
    res_norm = cblas_dnrm2(*la, residual, 1);
    resvec[k] = res_norm / rhs_norm;

    printf("Itération %d : Résidu normalisé = %.8f\n", k, resvec[k]);

    fprintf(file, "%d %.8f\n", k, resvec[k]);

    //Test de convergence.
    if (resvec[k] < *tol) 
    {
      *nbite = k + 1;
      break;
    }

    *nbite = k + 1;
  }

  fclose(file);
  free(MB);
  free(BX);
  free(residual);
}


double compute_error_Jacobi(double *X_numeric, double *X_exact, int la) 
{
  double *error = malloc(la * sizeof(double));
  for (int i = 0; i < la; i++) 
  {
    error[i] = X_numeric[i] - X_exact[i];
  }
  double error_norm = cblas_dnrm2(la, error, 1);
  free(error);
  return error_norm;
}

/*EXERCICE 9*/

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) 
{
  //Initialiser MB à zéro
  for (int i = 0; i < (*la) * (*kv); i++) 
  {
    MB[i] = 0.0;
  }

  //Copier la diagonale principale ainsi que la sous-diagonale dans MB
  for (int i = 0; i < *la; i++) 
  {
    //Copier la diagonale principale
    MB[i * (*kv) + *kl] = AB[i * (*lab) + *kl];

    //Copier la sous-diagonale (à la condition qu'elle existe)
    if (i > 0) 
    {
      MB[i * (*kv) + (*kl - 1)] = AB[i * (*lab) + (*kl - 1)];
    }
  }
}

void gauss_seidel_tridiag(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
  int kv = *kl + *ku + 1;     //Taille d'une ligne dans AB
  double *MB = malloc((*la) * kv * sizeof(double));   //Matrice triangulaire inférieure M
  double *residual = malloc(*la * sizeof(double));    //Vecteur des résidus
  double res_norm, rhs_norm;

  //Extraction de la matrice triangulaire inférieure M
  extract_MB_gauss_seidel_tridiag(AB, MB, lab, la, ku, kl, &kv);

  //Initialisation

  rhs_norm = cblas_dnrm2(*la, RHS, 1);
  if (rhs_norm == 0.0) {rhs_norm = 1.0;}

  FILE *file = fopen("résidus_gauss_seidel.txt", "w");
  if (file == NULL) 
  {
    perror("Erreur lors de l'ouverture du fichier");
    return;
  }

  for (int k = 0; k < *maxit; k++) 
  {
    //Itération de Gauss-Seidel
    for (int i = 0; i < *la; i++) 
    {
      //Calcul du résidu pour la ligne i
      double sigma = 0.0;

      //Contribution de la partie inférieure (j < i)
      if (i > 0) 
      {
        sigma += MB[i * kv + (*kl - 1)] * X[i - 1];
      }

      //Contribution de la partie supérieure (j > i)
      if (i < *la - 1) 
      {
        sigma += AB[i * (*lab) + (*kl + 1)] * X[i + 1];
      }

      //Mettre à jour X[i]
      X[i] = (RHS[i] - sigma) / MB[i * kv + *kl];
    }

    //Calcul du résidu

    cblas_dcopy(*la, RHS, 1, residual, 1); // residual = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, residual, 1);
    res_norm = cblas_dnrm2(*la, residual, 1);
    resvec[k] = res_norm / rhs_norm;

    printf("Itération %d : Résidu normalisé = %.8f\n", k, resvec[k]);

    fprintf(file, "%d %.8f\n", k, resvec[k]);

    //Test de convergence
    if (resvec[k] < *tol) 
    {
      *nbite = k + 1;
      break;
    }

    *nbite = k + 1;
  }
  fclose(file);

  free(MB);
  free(residual);
}

double calculate_error_Gauss_Seidel(double *X_numeric, double *X_exact, int la) 
{
  double *error = malloc(la * sizeof(double));
  for (int i = 0; i < la; i++) 
  {
    error[i] = X_numeric[i] - X_exact[i];
  }
  double error_norm = cblas_dnrm2(la, error, 1);
  free(error);
  return error_norm;
}



void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  //Initialisation des variables

  int kv = *kl + *ku + 1;   //Taille d'une ligne dans AB
  double *residual = malloc(*la * sizeof(double));    //Vecteur de résidus
  double *BX = malloc(*la * sizeof(double));    //Vecteur temporaire BX
  double rhs_norm, res_norm;
  double alpha_rich;    //Paramètre alpha de Richardson

  //Calcul de la norme du vecteur RHS (b)

  rhs_norm = cblas_dnrm2(*la, RHS, 1);
  if (rhs_norm == 0.0) { rhs_norm = 1.0; }

  *nbite = 0;  //Initialisation du compteur d'itérations

  //Idem qu'avant

  FILE *file = fopen("résidus_richardson_MB.txt", "w");
  if (file == NULL) 
  {
    perror("Erreur lors de l'ouverture du fichier");
    free(residual);
    free(BX);
    return;
  }

  for (int k = 0; k < *maxit; k++) 
  {
    //Calcul du résidu : r = b - Ax

    cblas_dcopy(*la, RHS, 1, BX, 1);    //BX = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, BX, 1);   //BX = RHS - AX
        
    //Calcul de alpha (facteur de Richardson)
    alpha_rich = richardson_alpha_opt(la);  //On peut utiliser la méthode définie pour alpha

    //Mise à jour de X : X = X + alpha * r
    cblas_daxpy(*la, alpha_rich, BX, 1, X, 1);

    //Calcul de la norme du résidu

    cblas_dcopy(*la, RHS, 1, residual, 1);  //residual = RHS
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, residual, 1);  //residual = RHS - AX
    res_norm = cblas_dnrm2(*la, residual, 1);  //Norme du résidu

    //Calcul du résidu normalisé

    resvec[k] = res_norm / rhs_norm;
    printf("Itération %d : Résidu normalisé = %.8f\n", k, resvec[k]);
    fprintf(file, "%d %.8f\n", k, resvec[k]);

    //Test de convergence

    if (resvec[k] < *tol) 
    {
      *nbite = k + 1;
      break;
    }
    *nbite = k + 1; 
  }

  //Fermer le fichier et libérer la mémoire

  fclose(file);
  free(residual);
  free(BX);
}

