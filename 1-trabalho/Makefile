#definição de variáveis
OBJ=main.o pdeSolver.o GaussSeidel.o Sor.o

all:	inicial pdeSolver clean

inicial: 
	@echo "Iniciando"

pdeSolver:$(OBJ)
	@gcc -o pdeSolver $(OBJ) -lm
main.o: pdeSolver.h main.c
	@gcc -c main.c
pdeSolver.o: pdeSolver.h pdeSolver.c
	@gcc -c pdeSolver.c
GaussSeidel.o: pdeSolver.h GaussSeidel.c
	@gcc -c GaussSeidel.c
Sor.o: pdeSolver.h Sor.c
	@gcc -c Sor.c

#remove arquivos objetos
clean:
	@echo "Finalizando"
	@rm -rf *.o 
