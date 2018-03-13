#include "pdeSolver.h"


double timestamp(void){
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}


void inicializa(double **X, double **B, int ny, int nx, double *fx, double *fy, double hx) {
    int i = 0,j;
    double superior, mul_PI;
 
 		mul_PI = 4.0*((3.1415)*(3.1415));
    superior =  sinh(2.0*3.1415);

    for(j = 1 ; j < (nx+1); ++j) {
      X[i][j] = sin(2.0*3.1415 * j * hx) * superior;
    }

    for (i = 1; i < (ny+1); ++i)											//tambem tinha 2 linhas a mais... mas ja fizemos a primeira... resto = 0
        for (j = 0 ; j < (nx+2); ++j) 							//nx colunas + 2 extras (a primeira e a ultima = 0)
          X[i][j] = 0.0;
    
    
    for( i=0; i < ny; ++i ) 
    		for(j = 0; j < nx; ++j)
    			B[i][j] = mul_PI * fy[i-1] * fx[j];

    return;
}


void impressao( double t_metodo, double t_residuo, double *residuo, char *metodo, FILE *saida, int maxiter, double hx, double hy, double **X, int nx, int ny) {
		int i, j;
    double y, x;

		fprintf(saida,"###########\n# Tempo Método %s: %lf\n", metodo, t_metodo );			//escrever o tempo medio
		fprintf(saida,"# Tempo Resíduo: %lf\n", t_residuo );
		fprintf(saida,"#\n# Norma do Residuo\n");

		for(i = 0; i < maxiter; ++i)				//escrever os valores do residuo de cada iteração
			fprintf(saida,"# i=%d: %.19lf\n", i, residuo[i]);
		fprintf(saida,"###########\n");
		
		//para o Gnuplot
		fprintf(saida, "set title \"Superficie da PDE para o metodo ");
		
  	if(strcmp(metodo,"GS") == 0)				//escrever o metodo
    	fprintf(saida, "GS\"\n");
  	else
    	fprintf(saida, "SOR\"\n");
    	
  	fprintf(saida, "set view 60,20\n");    //angulo de visão do gráfico
  	fprintf(saida, "set grid\n");
  	fprintf(saida, "splot '-' matrix nonuniform with lines\n");
  	fprintf(saida, "%d ", ny+1);
  	
  	for(x=0.0, j=0 ; j < (ny+2) ; x+=hx, ++j)						//valores de x
    	fprintf(saida, "%lf ", x);
  	fprintf(saida, "\n");
  	
																										  	// escreve a matriz e, no começo de cada linha, o respectivo valor y
  	for(i=0, y=1.0 ; i < (nx+2) ; ++i, y-=hy){
    	fprintf(saida, "%lf ", y);
   		for(j=0 ; j < (ny+2) ; ++j)
      	fprintf(saida, "%.19lf ", X[i][j]);
    	fprintf(saida, "\n");
  	}
  	fprintf(saida, "e\ne\n");
}
