Favor incluir neste diret�rio o arquivo contendo as posi��es iniciais e nome�-lo para in.dat

Como compilar
nvcc -ccbin clang-3.8 rendezvous.c -o rendezvous -Wno-deprecated-gpu-targets -lm -lOpenCL

Como executar
time optirun ./rendezvous N
