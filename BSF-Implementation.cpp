/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-M Skeleton
Module: BSF-Implementation.h (Implementation of BSF-M Skeleton)
Prefix: BI
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/
#include "BSF-Include.h"
#include "Problem-Parameters.h"		// Problem Parameters 
#include "Problem-Types.h"			// Problem Types 
#include "BSF-Types.h"				// BSF Types 
#include "BSF-Data.h"				// Problem Data 
#include "Problem-Forwards.h"		// Problem Parameters 
#include "BSF-Forwards.h"			// Function Forwards

using namespace std;

int main(int argc, char *argv[]) {
	BI_Init();

	char emptystring[] = "";
	char* message = emptystring;
	bool success = PI_bsf_Init(&BD_data, &message);
	if (!success) {
		if (BD_rank == BD_masterRank) {
			cout << "Error: PI_bsf_Init has returned False! ";
			if (strlen(message))
				cout << message;
			cout << endl;
		};
		MPI_Finalize();
		exit(1);
	};

	BI_MeasureTimeParameters();

	if (BD_rank == BD_masterRank)
		BI_Master();
	else
		BI_Worker();
	
	MPI_Finalize();
	return 0;
};
static void BI_Master() {// Master Process
	PI_bsf_ParametersOutput(BD_numOfWorkers, BD_data);
	BD_iterCount = 0;

	BD_t -= MPI_Wtime();
	do {
		BI_MasterMap(!BD_EXIT);
		BI_MasterReduce();
		BD_t_p -= MPI_Wtime();
		PI_bsf_ProcessResults(&BD_exit, &BD_data);
		BD_t_p += MPI_Wtime();
		BD_iterCount++;
		PI_bsf_IterOutput(BD_data, BD_iterCount, BD_t+MPI_Wtime());
	} while (!BD_exit);
	BD_t += MPI_Wtime();

	BI_MasterMap(BD_EXIT);

	BD_t_p /= (double)BD_iterCount;

	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		MPI_Irecv(&BD_time_w_Array[rank],sizeof(double),MPI_BYTE,rank,0,MPI_COMM_WORLD,&BD_request[rank]);
	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);

	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		BD_t_w += BD_time_w_Array[rank];
	BD_t_w /= (double)BD_iterCount;

	PI_bsf_ProblemOutput(BD_data, BD_iterCount,	BD_t, BD_t_L, BD_t_s, BD_t_r, BD_t_w, BD_t_p);
};

static void BI_Worker() {// Worker Process
	bool exit;

	while (true) {
		exit = BI_WorkerMap();
		if (exit) break;
		BI_WorkerReduce();
	};

	MPI_Send(&BD_t_w,sizeof(double),MPI_BYTE,BD_masterRank,0,MPI_COMM_WORLD);
};

static void BI_MasterMap(bool exit) {
#define APPENDIX(tail) (tail>0 ? 1 : 0)	
	int tail = BD_tailLength;
	int index = 0;

	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		PI_bsf_CopyData(&BD_data, &(BD_order[rank].data));
		BD_order[rank].exit = exit;
		BD_order[rank].index = index;
		BD_order[rank].length = BD_elemsPerWorker + APPENDIX(tail);
		MPI_Isend(
			&BD_order[rank],
			sizeof(BT_order_T),
			MPI_BYTE,
			rank,
			0,
			MPI_COMM_WORLD,
			&BD_request[rank]);
		index += BD_order[rank].length;
		tail--;
	};
	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);
};

static void BI_MasterReduce() {

	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		MPI_Irecv(
			&BD_reduceList[BD_order[rank].index],
			sizeof(PT_bsf_reduceElem_T)*BD_order[rank].length,
			MPI_BYTE,
			rank,
			0,
			MPI_COMM_WORLD,
			&BD_request[rank]);
	};

	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);

	PI_bsf_ProcessReduceList(BD_reduceList, &BD_data);

	//*debug*/cout << "BI_MasterReduce: approx = ";
	//*debug*/for (int j = 0; j < PP_N; j++) cout << BD_data.approximation[j] << "\t"; cout << endl;
};

static bool BI_WorkerMap() {
	MPI_Recv(
		&BD_order[BD_rank],
		sizeof(BT_order_T),
		MPI_BYTE,
		BD_masterRank,
		0,
		MPI_COMM_WORLD,
		&BD_status[BD_rank]);

	if (BD_order[BD_rank].exit)
		return BD_EXIT;

	BD_t_w -= MPI_Wtime();

	for (int i = BD_order[BD_rank].index; i < BD_order[BD_rank].index + BD_order[BD_rank].length; i++) 
		PI_bsf_MapF(&BD_mapList[i], &BD_reduceList[i], &BD_order[BD_rank].data);

	BD_t_w += MPI_Wtime();

	return !BD_EXIT;
};

static void BI_WorkerReduce() {
	MPI_Send(
		&BD_reduceList[BD_order[BD_rank].index],
		sizeof(PT_bsf_reduceElem_T)*BD_order[BD_rank].length,
		MPI_BYTE,
		BD_masterRank,
		0,
		MPI_COMM_WORLD);
};

static void BI_Init() {// Initialization
	int rc = MPI_Init(NULL, NULL);	// Starting MPI 
	if (rc != MPI_SUCCESS) {
		printf("Error starting MPI program. Terminating! \n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	};
	MPI_Comm_rank(MPI_COMM_WORLD, &BD_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &BD_size);

	if (BD_size < 2) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be > 1" << endl;
		MPI_Finalize();
		exit(1);
	};


	PI_bsf_AssignListSize(&BD_listSize);

	BD_mapList = (PT_bsf_mapElem_T*)malloc(BD_listSize * sizeof(PT_bsf_mapElem_T));
	PI_bsf_AssignMapList(BD_mapList, BD_listSize);

	BD_reduceList = (PT_bsf_reduceElem_T*)malloc(BD_listSize * sizeof(PT_bsf_reduceElem_T));

	if (BD_size > BD_listSize + 1) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be < Map List Size + 2 =" << BD_listSize + 2 << endl;
		MPI_Finalize();
		exit(1);
	};

	BD_masterRank = BD_size - 1;
	BD_numOfWorkers = BD_size - 1;
	BD_elemsPerWorker = BD_listSize / BD_numOfWorkers;
	assert(BD_elemsPerWorker > 0);
	BD_tailLength = BD_listSize - BD_elemsPerWorker*BD_numOfWorkers;
	BD_status = (MPI_Status*)malloc(BD_numOfWorkers * sizeof(MPI_Status));
	BD_request = (MPI_Request*)malloc(BD_numOfWorkers * sizeof(MPI_Request));
	BD_order = (BT_order_T*)malloc(BD_numOfWorkers * sizeof(BT_order_T));
	BD_time_w_Array = (double*)malloc(BD_numOfWorkers * sizeof(double));
};

static void BI_MeasureTimeParameters() {
	bool dummyByte;																
	BT_order_T dummyOrded;	
	PT_bsf_reduceElem_T *dummyReduceList;
	double startTime, stopTime;

	dummyReduceList = (PT_bsf_reduceElem_T*)malloc(BD_listSize * sizeof(PT_bsf_reduceElem_T));

	if (BD_rank == BD_masterRank) {

		MPI_Barrier(MPI_COMM_WORLD); // Measuring latency time 
		startTime = MPI_Wtime();
		MPI_Send(&dummyByte, sizeof(bool), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		stopTime = MPI_Wtime();
		BD_t_L = (stopTime - startTime);

		MPI_Barrier(MPI_COMM_WORLD); // Measuring send time/
		startTime = MPI_Wtime();																				
		for (int rank = 0; rank < BD_numOfWorkers; rank++) 														
			MPI_Isend(&dummyOrded, sizeof(BT_order_T), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);																									//
		MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);													
		stopTime = MPI_Wtime();																					
		BD_t_s = (stopTime - startTime) / BD_numOfWorkers - BD_t_L;

		MPI_Barrier(MPI_COMM_WORLD); // Measurement of receive time
		startTime = MPI_Wtime();																				
		for (int rank = 0; rank < BD_numOfWorkers; rank++)
			MPI_Irecv(dummyReduceList, sizeof(PT_bsf_reduceElem_T)*BD_listSize, MPI_BYTE, rank, 0, MPI_COMM_WORLD,
				&BD_request[rank]);
		MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);
		stopTime = MPI_Wtime();																					
		BD_t_r = (stopTime - startTime) / BD_numOfWorkers - BD_t_L;
	} else {
		MPI_Barrier(MPI_COMM_WORLD);// Measuring latency time
		if (BD_rank == 0)
			MPI_Recv(&dummyByte, sizeof(bool), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD, &BD_status[BD_rank]);

		MPI_Barrier(MPI_COMM_WORLD);// Measurement of send time 
		MPI_Recv(&dummyOrded, sizeof(BT_order_T), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD, &BD_status[BD_rank]);

		MPI_Barrier(MPI_COMM_WORLD);// Measurement of receive time 
		MPI_Send(dummyReduceList, sizeof(PT_bsf_reduceElem_T)*BD_listSize, MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
	};
};
