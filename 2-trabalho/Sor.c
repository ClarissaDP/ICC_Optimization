#include "pdeSolver.h"


void sor(int nx, int ny, int maxiter, FILE *saida) {
    int i, j, iter;
    double mulPI, grade, w;
    double hx, hy, hxD, hyD;
    double  X[(ny+2)*(nx+2)], B[ny*nx], residuo[maxiter];
    double time, t_residuo, t_metodo;


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
    /*for(i = 0; i < ny; ++i)
        fy[i] = sinh(2.0*3.1415 * (i+1) * hy);
  
    for(j = 0; j < nx ; ++j)
        fx[j] = sin(2.0 * 3.1415 * (j+1) * hx);*/
		
		/*Inicializar a matriz X (valores dos pontos = 0), bem como armazenar os valores de f(x,y) na matriz B*/
		//inicializa(X, B, ny, nx, hx, hy, residuo);
		
		/*
    for(iter = 0; iter < maxiter; ++iter) {
    		residuo[iter] = 0.0;
    		
        for(i = 1; i < (nx+1) ; ++i) {
            for(j = 1; j < (ny+1) ; ++j){
            	time = timestamp();
            		X[i*nx + j] = (1 - w)*X[i*nx + j];
                	X[i*nx + j] += w * ( grade *( (mulPI*fy[ny-i]*fx[j-1]) + hxD*(X[(i-1)*nx + j] + X[(i+1)*nx + j]) + hyD*(X[i*nx + (j-1)] + X[i*nx + (j+1)]) ));
                t_metodo += (timestamp () - time);
                
                time = timestamp();
                	residuo[iter] += fabs(( B[(i-1)*nx + (j-1)] - ( (-hxD)*(X[(i-1)*nx + j] + X[(i+1)*nx + j]) - hyD*(X[i*nx + (j-1)] + X[i*nx + (j+1)]) + ((2.0*hxD) + (2.0*hyD) + mulPI)*X[i*nx + j]) ));
                t_residuo += (timestamp () - time);
            }
        }
        residuo[iter] = sqrt(residuo[iter]);
    }
   */
    t_metodo = t_metodo/maxiter;
    t_residuo = t_residuo/maxiter;
	
    impressao(t_metodo, t_residuo, residuo, "SOR", saida, maxiter, hx, hy, X, nx, ny);		
		
    return;
}
