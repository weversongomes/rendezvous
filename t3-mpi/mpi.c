#include <mpi.h>
#include <stdio.h>
#include <math.h>

double x=0, y=0, z=0, xl0=0, yl0=0, zl0=0;
double w = 0.743737;
double ww;
int Tmax = 86400;

float brute_A (float y0, float xl0, float gama, float X, float vex, float vey, float gamavey_ww, float gama_w, float vex2_w)
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

float brute_B (float yl0, float X, float vey, float gamavex_ww, float gama_wpow, float vey_w) {
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

float brute_E (float y0, float xl0, float X, float vex) {
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

float dX (float t, float vex, float gama, float X, float A, float B, float E, float G, float vey2_w, float gama_wpow) {
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

        result2 += resultFn * exp(-(n * gama*t));
    }
    return result1 + result2;
}

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
    int NPI = 1;
    FILE *arq;
    char url[] = "in.dat";
    arq = fopen(url, "r");
    double var1;

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

    for(int np = 1; np <= NPI; np++) {

        if(arq == NULL) {
            printf("Erro, nao foi possivel abrir o arquivo\n");
            return 0;
        } else {
            fscanf(arq,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &var1, &var1, &var1, &x, &y, &z, &var1, &xl0, &yl0, &zl0, &var1, &var1, &var1, &var1, &$
        }

    for (int X=min; X<=max; X++) {
    printf("Printing from processor %s, rank %d out of %d processors, calculating for X %d\n", processor_name, world_rank, world_size, X);
        for(float Ve = 0.5f; Ve<=5.0f; Ve += 0.5f) {

            double vex, vey, vez;
            vex = vey = vez =Ve*Ve/3;
            double vex2_w = (2*vex)/w;
            double vey_w = vey/w;
            double vex3 = vex*3;
            double vey2_w = vex2_w;
            double vex4 = vex*4;

            for(int expGama = -14;  expGama<=2; expGama++) {
                double gama = pow(10, expGama);
                double gama_w = gama/w;
                double gamavex_ww = (gama*vex)/ww;
                double gamavey_ww = gamavex_ww;
                double gama_wpow = (gama/w)*(gama/w);
                double A = brute_A (y, xl0, gama, X, vex, vey, gamavey_ww, gama_w, vex2_w);
                double B = brute_B (yl0, X, vey, gamavex_ww, gama_wpow, vey_w);
                double E = brute_E (y, xl0, X, vex);
                double G = brute_G (x, yl0, X, vex, vey, vex3);

                for(int t = 0; t < Tmax; t++) {
                    dX(t, vex, gama, X, A, B, E, G, vey2_w, gama_wpow);
                    if (t == 0 && Ve == 0.5 && X == 1 && expGama == -14) {
                        printf("VALOR AMOSTRA: %lf\n", dX(t, vex, gama, X, A, B, E, G, vey2_w, gama_wpow));
                    }
                }
            }
        }
    }

    }

    // Finalize the MPI environment.
    MPI_Finalize();
}
