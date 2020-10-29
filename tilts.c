/*
	The TILTS algorithm implementation for HotSpot.
	Yongkui Han, yhan@ecs.umass.edu,
	Israel Koren, C. M. Krishna,
	Architecture and Real-Time Systems (ARTS) Laboratory, UMass at Amherst,
	Refer to "TILTS: A Lightweight Architectural-Level Thermal Simulation Method" for details.
	Or "Temptor: A Lightweight Runtime Temperature Monitoring Tool Using Performance Counters" at TACS2006.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "util.h"
#include "flp.h"
#include "tilts.h"

// TILTS = Time Invariant Linear Thermal System
// use the linear system theory to calculate the transient temperature fast.
double new_temp[MAX_NODES];
double temp2[MAX_NODES]; double temp3[MAX_NODES];


#if 1

// for the TILTS algorithm
//double * Aj, *Bx;
extern double Aj[24576]; double Bx[24576];

#endif

// Generate the matrices A and B for the TILTS algorithm
void generate_impulse_response_array(RC_model_t *model, flp_t * flp, double time_elapsed, double *A, double *B)
{
	int i,j; int n; int n_input, n_node;
	double power[1024], temp[1024];
	//char outfilename[128];
	//FILE * outfile = NULL;

	n = flp->n_units;
	n_input = n;
	n_node =NL* n + EXTRA;
	printf("n:  %d   n_node:  %d\n",flp->n_units,model->block->n_nodes);

	//generate pir array.
	for(i=0; i<n_input; i++) {
		
		//reset the power array and temp array.
		for(j=0; j<1024; j++) {
			power[j] = 0.0;
			temp[j] = AMBIENT_TEMP;
		}
		power[i] = 1.0;

		/* compute temperature	*/
		//compute_temp(power, temp, n, time_elapsed);
		compute_temp(model, power, temp, time_elapsed);

		for(j=0; j<n_node; j++) {
			//pir[i][j] = temp[j] - AMBIENT_TEMP;
			B[j*n_input+i] = temp[j] - AMBIENT_TEMP;
		}
	}

	// generate tir array.
	for(i=0; i<n_node; i++) {

		//reset the power array and temp array.
		for(j=0; j<1024; j++) {
			power[j] = 0.0;
			temp[j] = AMBIENT_TEMP;
		}
		temp[i] = AMBIENT_TEMP + 1.0;

		/* compute temperature	*/
		//compute_temp(power, temp, n, time_elapsed);
		compute_temp(model, power, temp, time_elapsed);

		for(j=0; j<n_node; j++) {
			//tir[i][j] = temp[j] - AMBIENT_TEMP;
			A[j*n_node+i] = temp[j] - AMBIENT_TEMP;
		}
	}

}

// For delta_t, we have A and B. This function calculates new A and B for 2*delta_t
void calculate_2timeAB(double *A, double *B, int n_node, int n_input)
{
	int i;
	double *A2, *B2;
	A2 = (double*)malloc(n_node * n_node * sizeof(double));
	B2 = (double*)malloc(n_node * n_input * sizeof(double));
	if(!A2 || !B2) {	//no enough memory!
		return;
	}	
	FILE *val;
	val=fopen("acc","w");
	// Calculate A2 = A * A
	matrix_mul(A, A, A2, n_node, n_node, n_node,val);
	printf("%f %f    ",A[0],A2[0]);
	// calculate A + I
	for(i=0; i<n_node; i++) {
		A[i*n_node+i] += 1.0;
	}
	
	// calculate B2 = ( A + I ) * B
	matrix_mul(A, B, B2, n_node, n_node, n_input,val);
	
	// Copy the new matrices back into A and B.
	memcpy(A, A2, n_node * n_node * sizeof(double));
	memcpy(B, B2, n_node * n_input * sizeof(double));
	
	if(A2) free(A2);
	if(B2) free(B2);
}

// For delta_t, we have A and B. 
// This function calculates new A and B for (2^n)*delta_t
void calculate_ntimeAB(int level, double *A, double *B, int n_node, int n_input)
{
	int i;
	for(i=0; i<level; i++) {
		calculate_2timeAB(A, B, n_node, n_input);
	}
}

// the input n_node, n_input is overwrited.
void allocate_AjBx_arrays(int n_node, int n_input)
{
	//n_node = N_NODE; n_input = N_INPUT;
	/*if(Aj == NULL) {
		Aj = (double*)malloc(n_node * n_node * sizeof(double));
	}
	if(Bx == NULL) {
		Bx = (double*)malloc(n_node * n_input * sizeof(double));	
	}*/
}

// Call this function before using the TILTS algorithm
void initialize_tilts(RC_model_t *model, flp_t * flp, double sampling_intvl, double *A, double *B)
{
	int n_input; int n_node;
	FILE *fp_A,*fp_B;
	fp_A=fopen("matrixA","w");
	fp_B=fopen("matrixB","w");
	//sampling_intvl = 0.000005;
	n_input = flp->n_units; n_node =NL* n_input + EXTRA;
	printf("initialization begins %f\n",sampling_intvl);
	if(sampling_intvl > 6e-6) {
		printf("initial 1\n");	
		// if the sampling interval is very large, accelerate the calculation of A and B
				double base_intvl = 1e-6;
				double times = sampling_intvl / base_intvl;
				int exp2 = 2; int num_exp2 = 1;
				while (exp2 < times) {
					exp2 *= 2; num_exp2 ++;
				}
        generate_impulse_response_array(model, flp, sampling_intvl/exp2, A, B);
        
        printf("matrix cal %d times \n",num_exp2);
        calculate_ntimeAB(num_exp2, A, B, n_node, n_input);
        
        
    }
    else {
    	printf("initial 2\n");
        generate_impulse_response_array(model, flp, sampling_intvl, A, B);
    }
	printf("A matrix is:\n"); 
	printf("B matrix is:\n"); 
	print_matrix(A,fp_A, n_node, n_node);
    print_matrix(B,fp_B, n_node, n_input);
	printf("initialization finished\n");
}

// print the values in the matrix

void print_matrix(double * a,FILE *fp, int m , int n)
{
  int i, j;
  for(i=0; i<m; i++) {
	for(j=0; j<n; j++) {
		fprintf(fp,"%.4f\t", a[i*n+j]);
	}
	fprintf(fp,"\n");
  }
}

// matrix multiply, c = a * b, a is MxN, b is NxL, c will be MxL
void matrix_mul(double * a, double *b, double *c, int m, int n, int l,FILE *fp)
{
	int i, j, k;
	double sum;
	sum=0;
	//FILE *fp;
	//fp=fopen("test","w");
	for(i=0; i<m; i++) {
		for(j=0; j<l; j++) {
			sum = 0;
			for(k=0; k<n; k++) {
				sum += a[i*n+k] * b[k*l+j];
				//fprintf(fp, "%f  %f  %f %d %d %d %f\n",sum,a[i*n+k],b[k*l+j],i,n,k,Bx[0]);	
			}
			c[i*l+j] = sum;
			//printf("%f  ",sum);
		}
		//printf("power %f  %d\n",Bx[0],i); 
	}
	//fclose(fp);
}

// perform Ax+Bu operation
void compute_temp_matrix(double *power, double *temp, double *A, double * B, int n_input, int n_node)
{
	int i;
	FILE *fp_A,*fp_B;
	//fp_A=fopen("test_A","w");
	//fp_B=fopen("test_B","w");
	//double temp2[N_NODE]; double temp3[N_NODE];
	//printf("power %f  \n",Bx[0]); 
	matrix_mul((double*)A, temp, temp2, n_node, n_node, 1,fp_A);
	//printf("power %f  \n",Bx[0]); 
	matrix_mul((double*)B, power, temp3, n_node, n_input, 1,fp_B);
	for(i=0; i<n_node; i++) {
		temp[i] = temp2[i] + temp3[i];
	}
	//printf("P: %f  T  %f  \n",temp3[0],temp2[0]);
}

// same as compute_temp() in HotSpot, but using the TILTS algorithm.
//void compute_temp_tilts(int n_input, int n_node, double *power, double *temp, double time_elapsed)
void compute_temp_tilts(double *temp, double *power, double time_elapsed,int n_input,int n_node, double *A, double *B)
{
	int j;
	//int n_input, n_node;
	double test_v=AMBIENT_TEMP;
	
	//n_input = model->block->n_units;
	//n_node = model->block->n_nodes;

	// first remove the ambient temperature.
	for(j=0; j<n_node; j++) {
		temp[j] -= AMBIENT_TEMP;
	}
	//printf("in tilts: %f\n",temp[0]);
	//printf("L2_left: %f  L2  %f   %f\n",temp[0],temp[1]);
	//for(j=0;j<n_node;j++){
		//printf("power %f  ",Bx[0]);
	//}

	compute_temp_matrix(power, temp, A, B, n_input, n_node);
	//printf("in tilts2: %f\n",temp[0]);
	//printf("after L2_left: %f  L2  %f\n",temp[0],temp[1]);
	// adding the ambient temperature back.
	for(j=0; j<n_node; j++) {
		temp[j] += AMBIENT_TEMP;
	}

}

