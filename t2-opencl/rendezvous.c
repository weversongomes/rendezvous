#include <stdio.h>
#include <stdlib.h>
     
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
     
#define MAX_SOURCE_SIZE (0x100000)
#define LENGTH (100)
#define X (100)
#define gama (17)
#define Ve (10)
    
float x=0, y=0, xl0=0, yl0=0;
     
int main(int argc, char *argv[]) {
    cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue command_queue = NULL;
    cl_mem d_a = NULL, d_b = NULL, d_c = NULL;
    cl_program program = NULL;
    cl_kernel kernel = NULL;
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
     
    FILE *fp;
    char fileName[] = "./rendezvous.cl";
    char *source_str;
    size_t source_size;
    float* h_a = (float*) calloc(LENGTH, sizeof(float));   // a vector
    float* h_b = (float*) calloc(LENGTH, sizeof(float));   // b vector
    float* h_c = (float*) calloc(LENGTH, sizeof(float));   // c vector
     
    /* Load the source code containing the kernel*/
    fp = fopen(fileName, "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);
    
    h_a[0] = 1;
    h_a[1] = 20;
    h_a[2] = 3;
    h_a[3] = 4;
    h_a[4] = 5;
    
    h_b[0] = 1;
    h_b[1] = 20;
    h_b[2] = 3;
    h_b[3] = 4;
    h_b[4] = 5;
    
    /* Get Platform and Device Info */
    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
     
    /* Create OpenCL context */
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
     
    /* Create Command Queue */
    command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
     
    /* Create Memory Buffer */
    d_a = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * LENGTH, NULL, &ret);
    d_b = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * LENGTH, NULL, &ret);
    d_c = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * LENGTH, NULL, &ret);
    // read/write = CL_MEM_READ_WRITE
     
    /* Create Kernel Program from the source */
    program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
    (const size_t *)&source_size, &ret);
     
    /* Build Kernel Program */
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    
    /* Create OpenCL Kernel */
    kernel = clCreateKernel(program, "main", &ret);
    
    //--------------------------------------
    FILE *arq;
    char url[] = "in.dat";
    arq = fopen(url, "r");
    float var1;
    int NPI; // numero de posicoes iniciais
	float done = 0; // usada para mostrar o percentual calculado
    
	if (argv[1] != NULL && atoi(argv[1]) >= 1) {
		NPI = atoi(argv[1]);
		printf("Calculando para %d posicoes iniciais.\n", NPI);
	} else {
		NPI = 1;
		printf("Calculando para 1 posicao inicial.\n");
	}

	printf("0%c concluido\r", 37);
	fflush(stdout);
	
    for(int np = 1; np <= NPI; np++) {
        
        if(arq == NULL) {
            printf("Erro, nao foi possivel abrir o arquivo\n");
            exit(EXIT_FAILURE);
        } else {
            fscanf(arq,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &var1, &var1, &var1, &x, &y, &var1, &var1, &xl0, &yl0, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1, &var1);
        }
    
        h_a[0] = x;
        h_a[1] = y;
        h_a[2] = xl0;
        h_a[3] = yl0;

        // Write a and b vectors into compute device memory
        ret = clEnqueueWriteBuffer(command_queue, d_a, CL_TRUE, 0, sizeof(float) * LENGTH, h_a, 0, NULL, NULL);
        if (ret != CL_SUCCESS) {
            printf("Error: 1! %d\n", ret);
            exit(1);
        }
        ret = clEnqueueWriteBuffer(command_queue, d_b, CL_TRUE, 0, sizeof(float) * LENGTH, h_b, 0, NULL, NULL);
        if (ret != CL_SUCCESS) {
            printf("Error: 2! %d\n", ret);
            exit(1);
        }
        
        /* Set OpenCL Kernel Parameters */
        /*ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&d_a);
        ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&d_b);
        ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&d_c);*/
        
        ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_a);
        ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_b);
        ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_c);
        
        /* Execute OpenCL Kernel */
        size_t global[3] = {10, 100, 17};
        //ret = clEnqueueTask(command_queue, kernel, 0, NULL,NULL);
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 3, NULL, global, NULL, 0, NULL, NULL);
        if (ret != CL_SUCCESS) {
            printf("Error: Failed to EnqueueNDRangeKernel! %d\n", ret);
            exit(1);
        }
        
        /* Copy results from the memory buffer */
        //ret = clEnqueueReadBuffer(command_queue, d_c, CL_TRUE, 0, MEM_SIZE * sizeof(char), h_c, 0, NULL, NULL);
        ret = clEnqueueReadBuffer(command_queue, d_c, CL_TRUE, 0, sizeof(float) * LENGTH, h_c, 0, NULL, NULL);
        if (ret != CL_SUCCESS) {
            printf("Error: Failed to read output array! %d\n", ret);
            exit(1);
        }
        
        /* Display Results */
        printf("VALOR AMOSTRA: %f\n", h_c[9]);
        /*for(int t = 0; t < 10; t++) {
            printf("%f\n", h_c[t]);
        }*/
        done = 100*np/NPI;
        if (done < 100) {
            printf("%.2f%c concluido\r", done, 37);
        } else {
            printf("%.2f%c concluido\n", done, 37);
        }
        fflush(stdout);
    }
    
    /* Finalization */
    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(d_a);
    ret = clReleaseMemObject(d_b);
    ret = clReleaseMemObject(d_c);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);
     
    free(source_str);
    free(h_a);
    free(h_b);
    free(h_c);
     
    return 0;
}
