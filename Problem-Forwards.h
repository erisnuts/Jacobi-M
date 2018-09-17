/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Jacobi method implementation using BSF-M Skeleton
Module: Problem-Forwards.h (BSF Function Forwards)
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/

#pragma once

//====================== BSF Forwards ===========================
void PI_bsf_AssignListSize(int* listSize);
void PI_bsf_AssignMapList(PT_bsf_mapElem_T* mapList, int listSize);
void PI_bsf_CopyData(PT_bsf_data_T* dataIn, PT_bsf_data_T* dataOut);
void PI_bsf_IterOutput(PT_bsf_data_T data, int iterCount, double elapsedTime);
bool PI_bsf_Init(PT_bsf_data_T* data, char** message);
void PI_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, PT_bsf_data_T* data);
void PI_bsf_ParametersOutput(int numOfWorkers, PT_bsf_data_T data);
void PI_bsf_ProblemOutput(PT_bsf_data_T data, int iterCount, double t, double t_L, double t_s,
	double t_r, double t_w, double t_p);
void PI_bsf_ProcessReduceList(PT_bsf_reduceElem_T *reduceList, PT_bsf_data_T *data);
void PI_bsf_ProcessResults(bool* exit, PT_bsf_data_T* data);

//====================== Problem Forwards ===========================
