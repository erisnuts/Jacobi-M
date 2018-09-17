/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Jacobi method implementation using BSF-M Skeleton
Module: Problem-Implementation.h (Implementation of the Problem)
Prefix: PI
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/

#include "Problem-Include.h"
#include "Problem-Types.h"			// Problem Types 
#include "Problem-Data.h"			// Problem Data 
#include "Problem-Parameters.h"		// Problem Parameters 
#include "Problem-Forwards.h"		// Function Forwards
using namespace std;
#define MIN(x,y) (x<y?x:y)

void PI_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, PT_bsf_data_T* data) {
	//*debug*/cout << "rowNo=" << mapElem->rowNo;
	reduceElem->coordinate = 0;
	for (int j = 0; j < PP_N; j++)
		reduceElem->coordinate += PD_Alpha[mapElem->rowNo][j] * data->approximation[j];
	reduceElem->coordinate += PD_beta[mapElem->rowNo];
	//*debug*/cout << " coordinate=" << reduceElem->coordinate << endl;
};

void PI_bsf_ProcessReduceList(PT_bsf_reduceElem_T *reduceList, PT_bsf_data_T *data) {
	for (int j = 0; j < PP_N; j++)
		data->approximation[j] = reduceList[j].coordinate;
};

void PI_bsf_ProcessResults(bool* exit, PT_bsf_data_T* data) {
	static int counter = 0;	// Iteration Counter
	counter++;
	if (counter < PP_ITER_COUNT)
		*exit = false;
	else
		*exit = true;
};

void PI_bsf_CopyData(PT_bsf_data_T* dataIn, PT_bsf_data_T* dataOut) {
	for (int j = 0; j < PP_N; j++)
		dataOut->approximation[j] = dataIn->approximation[j];
};

bool PI_bsf_Init(PT_bsf_data_T* data, char** message) {
	cout << setprecision(PP_PRECISION);

	for (int i = 0; i < PP_N; i++) { // Generating Matrix A
		for (int j = 0; j < i; j++)
			PD_A[i][j] = 1;
		PD_A[i][i] = (PP_FLOAT_POINT_TYPE) (i * 2);
		for (int j = i + 1; j < PP_N; j++)
			PD_A[i][j] = 0;
	};
	PD_A[0][0] = 1;/**/

/*	PD_A[0][0] = 5;		PD_A[0][1] = -2;	PD_A[0][2] = 3;
	PD_A[1][0] = -3;	PD_A[1][1] = 9;		PD_A[1][2] = 1;
	PD_A[2][0] = 2;		PD_A[2][1] = -1;	PD_A[2][2] = -7;/**/

	for (int i = 0; i < PP_N; i++) // Generating Vector of right parts
		PD_b[i] = (PP_FLOAT_POINT_TYPE) (i + 2 * i);
		PD_b[0] = 1;/**/

/*	PD_b[0] = -1;
	PD_b[1] = 2;
	PD_b[2] = 3;/**/

	for (int i = 0; i < PP_N; i++) { // Clculating reduced matrix Alpha
		for (int j = 0; j < PP_N; j++)
			PD_Alpha[i][j] = -PD_A[i][j] / PD_A[i][i];
		PD_Alpha[i][i] = 0;
	};

	for (int i = 0; i < PP_N; i++) // Clculating reduced vector beta
		PD_beta[i] = PD_b[i] / PD_A[i][i];

	for (int i = 0; i < PP_N; i++) // Generating coordinates of initial appriximation
		data->approximation[i] = PD_beta[i];
	return true;
};

void PI_bsf_AssignListSize(int* listSize) {
	*listSize = PP_N;
};

void PI_bsf_AssignMapList(PT_bsf_mapElem_T *mapList, int listSize) {
	for (int j = 0; j < listSize; j++) // First column has number 0
		mapList[j].rowNo = j;
};

void PI_bsf_ProblemOutput(PT_bsf_data_T data,	int iterCount, double t, double t_L, double t_s, 
	double t_r, double t_w, double t_p) {
	cout << "Iteration count = " << iterCount << endl;
	cout << "t_L = " << t_L << "\tt_s = " << t_s << "\tt_r = " << t_r << "\tt_w = " << t_w << "\tt_p = " << t_p << endl;
	cout << "Runtime = " << t << endl;
	cout << "Solution: ";
	for (int j = 0; j < MIN(PP_OUTPUT_LIMIT, PP_N); j++)
		cout << data.approximation[j] << "\t";
	cout << (PP_OUTPUT_LIMIT < PP_N? "..." : "") << endl << "___________________________________________________________________________" << endl;
};

void PI_bsf_IterOutput(PT_bsf_data_T data, int iterCount, double elapsedTime) {
/*	cout << iterCount << ": \t";
	for (int j = 0; j < PP_N; j++) 
		cout << data.approximation[j] << "\t\t";
	cout << endl;/**/
};

void PI_bsf_ParametersOutput(int numOfWorkers, PT_bsf_data_T data) {
	cout << "---------- Jacobi M -----------" << endl;
	cout << "Number of workers: " << numOfWorkers << endl;
	cout << "n = " << PP_N << endl;
/*	cout << "Initial approximation:\t";
	for (int i = 0; i < MIN(PP_OUTPUT_LIMIT, PP_N); i++)
		cout << data.approximation[i] << "\t";
	cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;

	cout << "---------- A -----------" << endl;
	for (int i = 0; i < MIN(PP_OUTPUT_LIMIT, PP_N); i++) {
		for (int j = 0; j < MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << PD_A[i][j] << "\t\t";
		cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	};
	cout << (PP_OUTPUT_LIMIT < PP_N ? "........." : "") << endl;

	cout << "---------- b -----------" << endl;
	for (int i = 0; i < MIN(PP_OUTPUT_LIMIT, PP_N); i++)
		cout << PD_b[i] << endl;
	cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;

	cout << "---------- Alpha -----------" << endl;
	for (int i = 0; i < MIN(PP_OUTPUT_LIMIT, PP_N); i++) {
		for (int j = 0; j < MIN(PP_OUTPUT_LIMIT, PP_N); j++)
			cout << PD_Alpha[i][j] << "\t\t";
		cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	};
	cout << (PP_OUTPUT_LIMIT < PP_N ? "........." : "") << endl;

	cout << "---------- beta -----------" << endl;
	for (int i = 0; i < MIN(PP_OUTPUT_LIMIT, PP_N); i++)
		cout << PD_beta[i] << endl;
	cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "___________________________________________________________________________" << endl;/**/
};