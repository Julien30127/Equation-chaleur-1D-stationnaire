#include "lib_poisson1D.h"

//Stockage de la matrice Poisson 1D

void poisson1D_CSR(int n, double *values, int *columns, int *row_ptr) 
{
    int num_nonzeros = 3 * (n - 1);     //Trois valeurs par ligne, sauf pour les bords (tridiagonalité)
    values = (double *)malloc(num_nonzeros * sizeof(double));
    columns = (int *)malloc(num_nonzeros * sizeof(int));
    row_ptr = (int *)malloc((n + 1) * sizeof(int));

    int idx = 0;
    for (int i = 0; i < n; i++) 
    {
        row_ptr[i] = idx;   // row_ptr[i] indique où commencent les éléments de la ligne i dans values et columns
        
        
        if (i > 0) 
        {
            //élément sous-diagonal

            values[idx] = -1;   
            columns[idx] = i - 1;
            idx++;
        }
        
        // élément diagonal

        values[idx] = 2;       
        columns[idx] = i;
        idx++;

        if (i < n - 1) 
        {
            //élément supra-diagonal

            values[idx] = -1;   
            columns[idx] = i + 1;
            idx++;
        }
    }
    row_ptr[n] = idx;  //Dernière ligne de row_ptr
}


void poisson1D_CSC(int n, double *values, int *rows, int *col_ptr) 
{
    int num_nonzeros = 3 * (n - 1);     //Trois valeurs par ligne, sauf pour les bords
    values = (double *)malloc(num_nonzeros * sizeof(double));
    rows = (int *)malloc(num_nonzeros * sizeof(int));
    col_ptr = (int *)malloc((n + 1) * sizeof(int));

    int idx = 0;
    for (int j = 0; j < n; j++) 
    {
        col_ptr[j] = idx;   //col_ptr[j] indique où commencent les éléments de la colonne j dans values et rows
        
        if (j > 0) 
        {
            //élément sous-diagonal

            values[idx] = -1;   
            rows[idx] = j - 1;
            idx++;
        }

        //élément diagonal

        values[idx] = 2;       
        rows[idx] = j;
        idx++;

        if (j < n - 1) 
        {
            //élément supra-diagonal

            values[idx] = -1;   
            rows[idx] = j + 1;
            idx++;
        }
    }
    col_ptr[n] = idx;   //Dernière colonne de col_ptr
}


//Produit Matrice Vecteur


void dcsrmv(int n, double *values, int *columns, int *row_ptr, double *x, double *y) 
{
    for (int i = 0; i < n; i++) 
    {
        //Initialisation du vecteur résultat
        y[i] = 0.0;  
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) 
        {
            y[i] += values[j] * x[columns[j]];
        }
    }
}


void dcscmv(int n, double *values, int *rows, int *col_ptr, double *x, double *y) 
{
    for (int j = 0; j < n; j++) 
    {
        //Initialisation du vecteur résultat
        y[j] = 0.0;  
        for (int i = col_ptr[j]; i < col_ptr[j + 1]; i++) 
        {
            y[rows[i]] += values[i] * x[j];
        }
    }
}


//À présent, réimplémentons les trois méthodes principales avec ce stockage


void richardson_sparse(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
    //Déclaration

    int kv = *kl + *ku + 1;
    double *residual = malloc(*la * sizeof(double));
    double *BX = malloc(*la * sizeof(double));
    double rhs_norm, res_norm, alpha;

    //Initialisation

    rhs_norm = cblas_dnrm2(*la, RHS, 1);
    if (rhs_norm == 0.0) rhs_norm = 1.0;
    *nbite = 0;

    //Calcul du produit matrice-vecteur pour Richardson avec le nouveau stockage
    //r = RHS - AB * X
    csr_matrix_vector_product(AB, RHS, X, residual, *la);

    FILE *file = fopen("résidus_richardson_sparse.txt", "w");
    if (file == NULL) 
    {
        perror("Erreur lors de l'ouverture du fichier");
        free(residual);
        free(BX);
        return;
    }

    for (int k = 0; k < *maxit; k++) {
        //Calcul de alpha (facteur de Richardson)
        alpha = 2.0 / (eigmin_poisson1D(*la) + eigmax_poisson1D(*la));  //optimal

        //Mise à jour de X
        cblas_daxpy(*la, alpha, residual, 1, X, 1);

        //Calcul de la norme du résidu avec le nouveau stockage

        csr_matrix_vector_product(AB, X, residual, *la);
        res_norm = cblas_dnrm2(*la, residual, 1);
        resvec[k] = res_norm / rhs_norm;

        //Test de convergence

        if (resvec[k] < *tol) 
        {
            *nbite = k + 1;
            break;
        }
        *nbite = k + 1;
    }

    fclose(file);
    free(residual);
    free(BX);
}



void jacobi_sparse(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{    
    double *residual = malloc(*la * sizeof(double));
    double rhs_norm, res_norm;

    rhs_norm = cblas_dnrm2(*la, RHS, 1);
    if (rhs_norm == 0.0) {rhs_norm = 1.0;}

    FILE *file = fopen("résidus_jacobi_sparse.txt", "w");
    if (file == NULL) {
        perror("Erreur lors de l'ouverture du fichier");
        free(residual);
        return;
    }

    for (int k = 0; k < *maxit; k++) 
    {
        //Calcul du nouveau X
        for (int i = 0; i < *la; i++) 
        {
            double sigma = 0.0;
            for (int j = *kl; j <= *ku; j++) 
            {
                if (i - j >= 0) sigma += AB[i * (*lab) + j] * X[i - j];
            }
            X[i] = (RHS[i] - sigma) / AB[i * (*lab) + *kl];
        }

        //Calcul du résidu avec le nouveau stockage

        csr_matrix_vector_product(AB, X, residual, *la);
        res_norm = cblas_dnrm2(*la, residual, 1);
        resvec[k] = res_norm / rhs_norm;

        //Test de convergence
        if (resvec[k] < *tol) {
            *nbite = k + 1;
            break;
        }
        *nbite = k + 1;
    }

    fclose(file);
    free(residual);
}

void gauss_seidel_sparse(double *AB, double *RHS, double *X, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite) 
{
    double *residual = malloc(*la * sizeof(double));
    double rhs_norm, res_norm;

    rhs_norm = cblas_dnrm2(*la, RHS, 1);
    if (rhs_norm == 0.0) {rhs_norm = 1.0;}

    FILE *file = fopen("résidus_gauss_seidel_sparse.txt", "w");
    if (file == NULL) {
        perror("Erreur lors de l'ouverture du fichier");
        free(residual);
        return;
    }

    for (int k = 0; k < *maxit; k++) 
    {
        //Itération de Gauss-Seidel
        for (int i = 0; i < *la; i++) {
            double sigma = 0.0;
            for (int j = *kl; j <= *ku; j++) 
            {
                if (i - j >= 0) sigma += AB[i * (*lab) + j] * X[i - j];
            }
            X[i] = (RHS[i] - sigma) / AB[i * (*lab) + *kl];
        }

        //Calcul du résidu avec le nouveau stockage

        csr_matrix_vector_product(AB, X, residual, *la);
        res_norm = cblas_dnrm2(*la, residual, 1);
        resvec[k] = res_norm / rhs_norm;

        //Test de convergence
        if (resvec[k] < *tol) 
        {
            *nbite = k + 1;
            break;
        }
        *nbite = k + 1;
    }

    fclose(file);
    free(residual);
}
