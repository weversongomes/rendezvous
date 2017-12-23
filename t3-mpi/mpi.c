#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double x=0, y=0, z=0, xl0=0, yl0=0, zl0=0;
/*double w = 0.743737;
double ww = 0.553144725169;
int Tmax = 86400;*/
double tStart, tEnd;

double brute_A (double y0, double xl0, double gama, double X, double vex, double vey, double gamavey_ww, double gama_w, double vex2_w)
{
    double w = 0.743737;
    double result = 0;
    double aux;
    double sum = 0;

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

double brute_B (double yl0, double X, double vey, double gamavex_ww, double gama_wpow, double vey_w) {
    double result = 0;
    double sum = 0;
    double aux;
    double w = 0.743737;
    
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

double brute_E (double y0, double xl0, double X, double vex) {
    double result = 0;
    double w = 0.743737;

    result -= 3*vex*log((X+1)/X);
    result += 6*w*y0 - 3*xl0;

    return result;
}

double brute_G (double x0, double yl0, double X, double vex, double vey) {
    double result = 0;
    double sum = 0;
    double aux;
    double w = 0.743737;

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

double dX (double t, double vex, double gama, double X, double A, double B, double E, double G, double vey2_w, double gama_wpow) {
    double w = 0.743737;
    //otimizacao
    double wt = w*t;

    double resultFn = 0.0f;
    double result1 = 2.0f * (A*sin(wt)-B*cos(wt))+E*t;

    double result2 = G;

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
    double w = 0.743737;
    double ww = 0.553144725169;
    int Tmax = 86400;

    int NPI; // numero de posicoes iniciais
    if (argv[1] != NULL && atoi(argv[1]) > 1) {
        NPI = atoi(argv[1]);
        printf("Calculando para %d posicoes iniciais.\n", NPI);
    } else {
        NPI = 1;
        printf("Calculando para 1 posicao inicial.\n");
    }

    int min = 1, max = 100;
    FILE *arq;
    char url[] = "in.dat";
    arq = fopen(url, "r");
    double var1;
    double done = 0;

  for(int np = 1; np <= NPI; np++) {

        if(arq == NULL) {
            printf("Erro, nao foi possivel abrir o arquivo\n");
            return 0;
        } else {
            fscanf(arq,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &var1, &var1, &var1, &x, &y, &z, &var1, &xl0, &yl0, &zl0, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1);
        }
        
    for (int X=min; X<=max; X++) {
        printf("X = %d\n", X);
        for(double Ve = 0.5f; Ve<=5.0f; Ve += 0.5f) {

            double vex, vey, vez;
            vex = vey = vez =Ve*Ve/3;
            double vex2_w = (2*vex)/w;
            double vey_w = vey/w;
            double vex3 = vex*3;
            double vey2_w = vex2_w;
            double vex4 = vex*4;
            
            double E = brute_E (y, xl0, X, vex);
            double G = brute_G (x, yl0, X, vex, vey);

            for(int expGama = -14;  expGama<=2; expGama++) {
                double gama = pow(10, expGama);
                double gama_w = gama/w;
                double gamavex_ww = (gama*vex)/ww;
                double gamavey_ww = gamavex_ww;
                
                double gama_wpow = (gama/w)*(gama/w);
                double A = brute_A (y, xl0, gama, X, vex, vey, gamavey_ww, gama_w, vex2_w);
                double B = brute_B (yl0, X, vey, gamavex_ww, gama_wpow, vey_w);

                for(int t = 0; t < Tmax; t++) {
                    dX(t, vex, gama, X, A, B, E, G, vey2_w, gama_wpow);
                    
                    if (t == 0 && Ve == 0.5 && X == 1 && expGama == -14) {
                        double result = dX(t, vex, gama, X, A, B, E, G, vey2_w, gama_wpow);
                        printf("VALOR AMOSTRA: %lf\n", result);
                    }
                }
            }
        }
        //fflush(stdout);
    }
  }
}
