#include "pdeSolver.h"


void sor(int nx, int ny, int maxiter, FILE *saida) {
    int i, j, iter;
    double fx[nx], fy[ny], mulPI, grade, w;
    double hx, hy, hxD, hyD;
    double **X, **B, residuo[maxiter];
    double time, t_residuo, t_metodo;


   
    /* Alocações dinâmicas */
    if(!( X = (double**)malloc((ny+2)*sizeof(double)) )) {
        fprintf(stderr, "Erro na alocação da matriz X.\n");
        return;
    }   
    for(i = 0 ; i < (ny+2) ; ++i)
        if(!( X[i] = (double*)malloc((nx+2)*sizeof(double)) )) {
          fprintf(stderr, "Erro na alocação da matriz X.\n");
            return;
        }   


    if(!( B = (double**)malloc( nx *sizeof(double)) )) {
        fprintf(stderr, "Erro na alocação da matriz B.\n");
        return;
    }   
    for(i = 0 ; i < nx ; ++i)
        if(!( B[i] = (double*)malloc( ny *sizeof(double)) )) {
          fprintf(stderr, "Erro na alocação da matriz B.\n");
          return;
        }   



		/*Aqui estão variaveis/valores importantes para o metodo*/
		hx = (2.0/(nx+1));
		hy = (1.0/(ny+1));
		t_metodo = 0.0;
		t_residuo = 0.0;
		w = 2.0 - (hx+hy)/2.0;                                  //fator de relaxação

        hxD = (hx*hx);                                    //hxD = hx²
        hyD = (hy*hy);                                    //hyD = hy²

        mulPI = 4.0*((3.1415)*(3.1415));                          //4π²
        grade = (1.0 / ((2.0/hxD) + (2.0/hyD) + mulPI));          // 1/(2/hxD + 2/hyD + 4π² )

		    hxD = 1.0/hxD;                                    //hxD = 1/hx²
		    hyD = 1.0/hyD;                                    //hyD = 1/hy²


	/* Salvar todos os valores da função f(x,y) para todos os x e todos os y */
    for(i = 0; i < ny; ++i)
        fy[i] = sinh(2.0*3.1415 * (i+1) * hy);
  
    for(j = 0; j < nx ; ++j)
        fx[j] = sin(2.0 * 3.1415 * (j+1) * hx);
		
		/*Inicializar a matriz X (valores dos pontos = 0), bem como armazenar os valores de f(x,y) na matriz B*/
		inicializa(X, B, ny, nx, fx, fy, hx);
		
		
    for(iter = 0; iter < maxiter; ++iter) {
    		residuo[iter] = 0.0;
    		
        /* Sor na fórmula */
        for(i = 1; i < (nx+1) ; ++i) {
            for(j = 1; j < (ny+1) ; ++j){
            		/*metodo*/
            	time = timestamp();
            		X[i][j] = (1 - w)*X[i][j];
                	X[i][j] += w * ( grade *( (mulPI*fy[ny-i]*fx[j-1]) + hxD*(X[i-1][j] + X[i+1][j]) + hyD*(X[i][j-1] + X[i][j+1])) );
                t_metodo += (timestamp () - time);
                
                /*residuo*/
                time = timestamp();
                	residuo[iter] += fabs(( B[i-1][j-1] - ( (-hxD)*(X[i-1][j] + X[i+1][j]) - hyD*(X[i][j-1] + X[i][j+1]) + ((2.0*hxD) + (2.0*hyD) + mulPI)*X[i][j]) ));
                t_residuo += (timestamp () - time);
            }
        }
        residuo[iter] = sqrt(residuo[iter]);
    }
   
    t_metodo = t_metodo/maxiter;
    t_residuo = t_residuo/maxiter;
	
    impressao(t_metodo, t_residuo, residuo, "SOR", saida, maxiter, hx, hy, X, nx, ny);		
		
    return;
}
