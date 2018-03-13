#include "pdeSolver.h"


int main(int argc, char *argv[]) {
  int maxIter, i;
  int nx, ny;
  char m[3];		
  FILE *saida;				

	if( argc < 11 ) {
		printf("ERRO.Uso correto:\n%s -nx <Nx> -ny <Ny> -i <maxIter> -m <gs | sor> -o arquivo_saida\n", argv[0]);
		return (-1);
	}

	for (i = 1; i < 11 ; ++i) {
		if ( strcmp(argv[i],"-nx") == 0 )
			nx = atof(argv[++i]);
		else if ( strcmp(argv[i],"-ny") == 0 )
			ny = atof(argv[++i]);
		else if ( strcmp(argv[i],"-i") == 0 )
			maxIter = atof(argv[++i]);
		else if ( strcmp(argv[i],"-m") == 0 )
			strcpy(m,argv[++i]);
		else if ( strcmp(argv[i],"-o") == 0){
			saida = fopen(argv[++i],"w");
			if (saida == NULL){
				fprintf(stderr, "Erro ao abrir o arquivo de saída.\n");
				return (-1);
			}	
		}
		else
			puts("Argumento Inválido, melhor checar isso pessoa!");
	}


	/* Caso método escolhido seja gs -> Gauss Seidel */
	if (strcmp(m,"gs") == 0) 
		gaussSeidel(nx,ny,maxIter,saida);
	

	/* Caso método escolhido seja sor -> Successive Over Relaxation */
	else
		sor(nx,ny,maxIter,saida);
	

	fclose(saida);

  return 0;
}
