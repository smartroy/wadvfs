/*
	The TILTS algorithm implementation for HotSpot.
	Yongkui Han, yhan@ecs.umass.edu,
	Israel Koren, C. M. Krishna,
	Architecture and Real-Time Systems (ARTS) Laboratory, UMass at Amherst,
	Refer to "TILTS: A Lightweight Architectural-Level Thermal Simulation Method" for details.
	Or "Temptor: A Lightweight Runtime Temperature Monitoring Tool Using Performance Counters" at TACS2006.
 */


#ifndef __TILTS_H
#define __TILTS_H

#include <stdio.h>
//#include "RC.h"
#include "temperature.h"
#include "temperature_grid.h"
#include "temperature_block.h"

// yongkui, the following is devoted to the Linear system based fast transient thermal simulation.

//#define AMBIENT_TEMP	(273.15 + 40.0)
//#define AMBIENT_TEMP	(273.15 + 45.0)
#define AMBIENT_TEMP (273.15 + T_INIT)
// FOR the alpha processor
#define N_INPUT		18
//#define N_NODE		64			// for original ev6.flp without rim
//#define N_NODE		76		// for ev6.flp with 4 RIM blocks
#define N_NODE		97		// for ev6.flp with 11 RIM blocks

//for the pentium pro processor
//#define N_INPUT		16
//#define N_NODE		58

#define MAX_NODES		128
#define MAX_LEVELS		16

void generate_impulse_response_array(RC_model_t *model, flp_t * flp, double time_elapsed, double *A, double *B);
void calculate_2timeAB(double *A, double *B, int n_node, int n_input);
void calculate_ntimeAB(int level, double *A, double *B, int n_node, int n_input);

void allocate_AjBx_arrays(int n_node, int n_input);

void initialize_tilts(RC_model_t *model, flp_t * flp, double sampling_intvl,double *A, double *B);

//void print_matrix(double * a, FILE *fp,int m , int n);
void matrix_mul(double * a, double *b, double *c, int m, int n, int l,FILE *fp);

void compute_temp_matrix(double *power, double *temp, double *A, double * B, int n_input, int n_node);

//void compute_temp_tilts(int n_input, int n_node, double *power, double *temp, double time_elapsed);
void compute_temp_tilts(double *temp, double *power, double time_elapsed,int n_input,int n_node,double *A, double *B);
void print_matrix(double * a,FILE *fp, int m , int n);
#endif

