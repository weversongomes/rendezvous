Favor incluir neste diretório o arquivo contendo as posições iniciais e nomeá-lo para in.dat

Como compilar
nvcc -ccbin clang-3.8 rendezvous.c -o rendezvous -Wno-deprecated-gpu-targets -lm -lOpenCL

Como executar
time optirun ./rendezvous N
