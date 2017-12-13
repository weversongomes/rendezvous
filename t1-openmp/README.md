# Computacao-de-alto-desempenho

ATENCAO
Favor incluir o arquivo contendo as posicoes iniciais neste mesmo diretorio e nomea-lo para in.dat

Como compilar
gcc -Wall -O3 x.c -o x -lm -fopenmp

Como executar
time ./x N n
lista de parametros
	N - Numero de linhas para serem lidas do arquivo de entrada in.dat
	n - Numero de nucleos que serao utilizados
