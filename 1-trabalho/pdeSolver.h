#ifndef _pdeSolver_
#define _pdeSolver_

#include <stdio.h>                                                         
#include <stdlib.h>                                                        
#include <string.h> 
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

#endif


double timestamp(void);
void impressao( double, double, double *, char *, FILE *, int, double, double, double**, int, int);
void inicializa(double**, double**, int, int, double *, double *, double);


/*------------------------------------------------------------------------------*/
/* Declarações de funções do Método de Gauss Seidel															*/
/*------------------------------------------------------------------------------*/

void gaussSeidel(int, int, int, FILE *);


/*------------------------------------------------------------------------------*/
/* Declarações de funções do Método Successive Over Relaxation									*/
/*------------------------------------------------------------------------------*/
void sor(int, int, int, FILE *);

