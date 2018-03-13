#definição de variáveis
OBJ=main.o pdeSolver.o GaussSeidel.o Sor.o

all:	inicial pdeSolver clean

inicial: 
	@echo "Iniciando"

pdeSolver:$(OBJ)
	@gcc -pthread -o pdeSolver $(OBJ) -DLIKWID_PERFMON -L/home/soft/likwid/lib -llikwid -lm
	@#@gcc -pthread -o pdeSolver $(OBJ) -DLIKWID_PERFMON -L/usr/local/lib -llikwid -lm
main.o: pdeSolver.h main.c
	@gcc -c main.c
pdeSolver.o: pdeSolver.h pdeSolver.c
	@gcc -c pdeSolver.c
GaussSeidel.o: pdeSolver.h GaussSeidel.c
	@gcc -c GaussSeidel.c #-pthread -DLIKWID_PERFMON -L/usr/local/lib -llikwid -lm
Sor.o: pdeSolver.h Sor.c
	@gcc -c Sor.c

#remove arquivos objetos
clean:
	@echo "Finalizando"
	@rm -rf *.o 
