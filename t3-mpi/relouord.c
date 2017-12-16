#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    
    int min, max;
    if (world_rank == 0) {
        min = 1;
        max = 25;
    } else if (world_rank == 1) {
        min = 26;
        max = 50;
    } else if (world_rank == 2) {
        min = 51;
        max = 75;
    } else if (world_rank == 3) {
        min = 76;
        max = 100;
    }
    
    for (int X=min; X<=max; X++) {
    // Print off a hello world message
    printf("Printing from processor %s, rank %d out of %d processors, calculating for X %d\n", processor_name, world_rank, world_size, X);
        for(float Ve = 0.5f; Ve<=5.0f; Ve += 0.5f) {
            for(int expGama = -14;  expGama<=2; expGama++) {
                
            }
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}

float brute_A (float y0, float xl0, float gama, float X, float vex, float vey, float gamavey_ww, float gama_w, float vex2_w, float w)
{
    float result = 0;
    float aux;
    float sum = 0;

    result = (2*xl0)/w - 3*y0 + ((2*vex)/w) * log((X+1)/X);

    // Calculo do somatorio
    for (int n = 1; n <= 20; n++) {
        aux = (1/(n*pow(X, n)))*(1/(1+(n*gama_w)*(n*gama_w)))*((vex2_w)+(n*gamavey_ww));
        if (n%2 == 0) {//iteracao par
            aux = -aux;
        }
        sum += aux;
    }
    result -= sum;

    return result;
}

float brute_B (float yl0, float X, float vey, float gamavex_ww, float gama_wpow, float vey_w, float w) {
    float result = 0;
    float sum = 0;
    float aux;

    result = yl0/w + (vey/w)*log((X+1)/X);

    // Calculo do somatorio
    for (int n = 1; n <= 20; n++) {
        aux = (1/(n*pow(X,n)))*(1/(1+(n*n*gama_wpow)))*(vey_w + (n*gamavex_ww));
        if (n%2 == 0) {//iteracao par
            aux = -aux;
        }
        sum += aux;
    }

    result+= sum;
    
    return result;
}

float brute_E (float y0, float xl0, float X, float vex, float w) {
    float result = 0;

    result -= 3*vex*log((X+1)/X);
    result += 6*w*y0 - 3*xl0;

    return result;
}

float brute_G (float x0, float yl0, float X, float vex, float vey, float w) {
    float result = 0;
    float sum = 0;
    float aux;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;
    
    for (int n = 1; n <= 20; n++) {
        aux = (vex*3)/(n*n*pow(X,n)*w);

        if (n%2 == 0) {
            aux = -aux;
        }
        sum +=aux;
    }

    result-=sum;

    return result;
}

float dX (float t, float vex, float gama, float X, float A, float B, float E, float G, float vey2_w, float gama_wpow, float w) {
    //otimizacao
    float wt = w*t;

    float resultFn = 0.0f;
    float result1 = 2.0f * (A*sin(wt)-B*cos(wt))+E*t;
    
    float result2 = G;

    for (int n = 1; n <= 20; n++) {
        // brute_F
        resultFn = (1/(n*pow(X,n)))*(vey2_w + (vex*4)/(n*gama))/(1+(n*n*gama_wpow));

        if (n%2 == 0) {
            resultFn = -resultFn;
        }
        resultFn -= vex/(n*gama);
        //brute_F

        result2 += resultFn * pow(M_E, -(n * gama*t));
    }
    return result1 + result2;
}
