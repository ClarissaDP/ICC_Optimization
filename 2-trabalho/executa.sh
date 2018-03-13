#!/bin/bash

echo "tamanhoMatriz caminhoDiretórioSaida"

# Gera 10 execuções (para FLOPS_DP, CACHE e MEM cada) para um determinado tamanho inserido e coloca no caminho indicado.
# Ver depois qual foi a melhor execução para cada.

tam=$1
caminho=$2


for ((i = 0; i < 10; i++)) do
  
  likwid-perfctr -C 6 -g FLOPS_DP -m ./pdeSolver -nx $tam -ny $tam -i 10 -m gs -o "$caminho"saida-"$tam"-"$i".txt > "$caminho""$tam"-"$i"-FLOPS_DP.txt
  likwid-perfctr -C 6 -g CACHE -m ./pdeSolver -nx $tam -ny $tam -i 10 -m gs -o "$caminho"saida-"$tam"-"$i".txt > "$caminho""$tam"-"$i"-CACHE.txt
  likwid-perfctr -C 6 -g MEM -m ./pdeSolver -nx $tam -ny $tam -i 10 -m gs -o "$caminho"saida-"$tam"-"$i".txt > "$caminho""$tam"-"$i"-MEM.txt

done

