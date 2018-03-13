#include "pdeSolver.h"


void gaussSeidel(int nx, int ny, int maxiter, FILE *saida) {
    int i, j, iter;
    double fx[nx], fy[ny], mulPI, grade, aux1, aux2, aux3, aux4;
    double hx, hy, hxD, hyD;
    double *X, *B, residuo[maxiter];
    double time, t_residuo, t_metodo;

		X = (double *) malloc ((ny+2)*(nx+2)*sizeof(double));
		B = (double *) malloc (ny*nx*sizeof(double));

		/* Aqui estão variaveis/valores importantes para o metodo */
		hx = (2.0/(nx+1));
		hy = (1.0/(ny+1));
		t_metodo = 0.0;
		t_residuo = 0.0;

        hxD = (hx*hx);                                        //hxD = hx²
        hyD = (hy*hy);                                        //hyD = hy²

        mulPI = 4.0*((3.1415)*(3.1415));                              //4π²
        grade = (1.0 / ((2.0/hxD) + (2.0/hyD) + mulPI));              // 1/(2/hxD + 2/hyD + 4π² )
        
        hxD = 1.0/hxD;                                        //hxD = 1/hx²
        hyD = 1.0/hyD;                                        //hyD = 1/hy²
	

	/* Salvar todos os valores da função f(x,y) para todos os x e todos os y */
    for(i = 0; i < ny; ++i)
        fy[i] = sinh(2.0*3.1415 * (i+1) * hy);
  
    for(j = 0; j < nx ; ++j)
        fx[j] = sin(2.0 * 3.1415 * (j+1) * hx);
		
    aux2 = ((2.0*hxD) + (2.0*hyD) + mulPI);  
		/* Inicializar a matriz X (valores dos pontos = 0), bem como armazenar os valores de f(x,y) na matriz B*/
		inicializa(X, B, ny, nx, fx, fy, hx, hy, residuo);
		
		likwid_markerInit();
    likwid_markerStartRegion("Método");
    aux3 = X[nx + 2];
    for(iter = 0; iter < maxiter; ++iter) {
    		residuo[iter] = 0.0;
    		
        //Gauss na fórmula /
        for(i = 1; i <= ny ; ++i) {

            /* Cobrir não múltiplos */
            for(j = 1; j <= (nx % 4); ++j) {
            
            		aux1 = hxD*(X[(i-1) * nx +j] + X[(i+1)*nx + j]) + hyD*(aux3 + X[i*nx +j+1]);
            		aux3 = grade *( B[(i-1) * nx + j - 1] + aux1);
            		X[i * nx + j] = aux3;
            
            }

            for(; j <= nx ; j += 4) {

              /*metodo*/
            		time = timestamp();
            		
            		aux1 = hxD*(X[(i-1) * nx +j] + X[(i+1)*nx + j]) + hyD*(aux3 + X[i*nx +j+1]);
            		aux3 = grade *( B[(i-1) * nx + j - 1] + aux1);
            		X[i * nx + j] = aux3;
            		
            		aux1 = hxD*(X[(i-1) * nx +j + 1] + X[(i+1)*nx + j + 1]) + hyD*(aux3 + X[i*nx +j+2]);
            		aux3 = grade *( B[(i-1) * nx + j] + aux1);
            		X[i * nx + j + 1] = aux3;
            		
            		aux1 = hxD*(X[(i-1) * nx +j + 2] + X[(i+1)*nx + j + 2]) + hyD*(aux3 + X[i*nx +j+3]);
            		aux3 = grade *( B[(i-1) * nx + j + 1] + aux1);
            		X[i * nx + j + 2] = aux3;
            		
            		aux1 = hxD*(X[(i-1) * nx +j + 3] + X[(i+1)*nx + j + 3]) + hyD*(aux3 + X[i*nx +j+4]);
            		aux3 = grade *( B[(i-1) * nx + j + 2] + aux1);
            		X[i * nx + j + 3] = aux3;
            			
                t_metodo += (timestamp () - time);  
            }
         //residuo
         for(i = 1; i <= ny ; ++i) {
         		aux3 = X[i*(nx +2)];							//é o anterior
         		aux1 = X[i*(nx + 2) + 1];					//o atual
            aux4 = X[i*(nx + 2) + j + 2];			//próximo
            
            /* Cobrir não múltiplos */
            for(j = 1; j <= (nx % 4); ++j) {
                
                	residuo[iter] += fabs(  - ( aux2 * aux1 - (hxD*(X[(i-1) * nx +j] + X[(i+1)*nx + j]) + hyD*(aux3 + aux4))));
                	aux3 = aux1;
                	aux1 = aux4;
                	aux4 = X[i*nx + j + 2];
            }
            
            for(; j <= nx ; j += 4) {
            
              /*residuo*/
                time = timestamp();
                
                	residuo[iter] += fabs(  - ( aux2 * aux1 - (hxD*(X[(i-1) * nx +j] + X[(i+1)*nx + j]) + hyD*(aux3 + aux4))));
                	aux3 = aux1;
                	aux1 = aux4;
                	aux4 = X[i*nx + j + 2];
                	
                	residuo[iter] += fabs(  - ( aux2 * aux1 - (hxD*(X[(i-1) * nx +j + 1] + X[(i+1)*nx + j + 1]) + hyD*(aux3 + aux4))));
                	aux3 = aux1;
                	aux1 = aux4;
                	aux4 = X[i*nx + j + 2];
                	
                	residuo[iter] += fabs(  - ( aux2 * aux1 - (hxD*(X[(i-1) * nx +j + 2] + X[(i+1)*nx + j + 2]) + hyD*(aux3 + aux4))));
                	aux3 = aux1;
                	aux1 = aux4;
                	aux4 = X[i*nx + j + 2];
                	
                	residuo[iter] += fabs(  - ( aux2 * aux1 - (hxD*(X[(i-1) * nx +j + 3] + X[(i+1)*nx + j + 3]) + hyD*(aux3 + aux4))));
                	aux3 = aux1;
                	aux1 = aux4;
                	aux4 = X[i*nx + j + 2];
                	
                t_residuo += (timestamp () - time);
            }
         }
        }
        residuo[iter] = sqrt(residuo[iter]);
   	}
   	
   	
   	likwid_markerStopRegion("Método");
    likwid_markerClose();
   
   
   
    t_metodo = t_metodo/maxiter;
    t_residuo = t_residuo/maxiter;
	
    impressao(t_metodo, t_residuo, residuo, "GS", saida, maxiter, hx, hy, X, nx, ny);			

    free(X);
    free(B);

    return;
}
