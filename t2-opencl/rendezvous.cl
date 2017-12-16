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
    
__kernel void main(__global float* a, __global float* b, __global float* c)
{
    int Ve = get_global_id(0);
    int X = get_global_id(1);
    int expGama = get_global_id(2);
    float w = 0.743737;
    float gama;
    float vex, vey;
    float vex2_w, vey2_w;
    //otimizacao -------------------------------
    //brute_A
    float gama_w;
    //brute_B
    float gamavex_ww;
    //dy vy A
    float gamavey_ww;
    //vx vy vz dx dy dz B H I
    float gama_wpow;
    float A = 0.0f, B = 0.0f, E = 0.0f, G = 0.0f;
    float VeReal = 0.0f;
    float response = 0.0f;
    //------------------------------------------
    c[0] = get_global_size(0);
    c[1] = get_global_size(1);
    c[2] = get_global_size(2);
    c[3] = a[0];
    c[4] = a[1];
    c[5] = a[2];
    c[6] = a[3];
    if(Ve < 10 && X < 100 && expGama < 17) {
        VeReal = 0.5f * (Ve + 1);
        gama = pown(10.0f, expGama - 14);
        vex = vey = VeReal*VeReal/3;
        vex2_w = vey2_w = (2*vex)/w;
        gama_w = gama/w;
        gamavex_ww = gamavey_ww = (gama*vex)/(w*w);
        gama_wpow = (gama/w)*(gama/w);
        
        A = brute_A (a[1], a[2], gama, X+1, vex, vey, gamavey_ww, gama_w, vex2_w, w);
        B = brute_B (a[3], X+1, vey, gamavex_ww, gama_wpow, vex2_w, w);
        E = brute_E (a[1], a[2], X+1, vex, w);
        G = brute_G (a[0], a[3], X+1, vex, vey, w);
        for(int t = 0; t < 86400; t++) {
            if (t == 0 && Ve == 0 && X == 0 && expGama == 16) {
                c[9] = dX(t, vex, gama, X+1, A, B, E, G, vey2_w, gama_wpow, w);
                printf("Printing from GPU\n");
            } else {
                dX(t, vex, gama, X+1, A, B, E, G, vey2_w, gama_wpow, w);
            }
        }
        
        if(Ve == 5 && X == 56 && expGama == 16) {
            c[7] = 0.0f;
        } else if(Ve == 0 && X == 0 && expGama == 0) {
            c[8] = A;
        }
    }
}
