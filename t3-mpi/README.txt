mpicc relouord.c -o relouord -std=c99 -fopenmp
mpirun -np 4 -machinefile /home/weverson/machinefile.txt /home/weverson/relouord N M
    N é o número de posições iniciais
    M é o número de núcleos usados por processador
