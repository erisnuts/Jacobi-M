/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Jacobi method implementation using BSF-M Skeleton
Module: Problem-Types.h (BSF Types)
Prefix: PT
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/

#pragma once						
#include "Problem-Parameters.h"		// Problem Parameters 

//=========================== Problem Types =========================
typedef PP_FLOAT_POINT_TYPE PT_point_T[PP_N];	// Point in n-Dimensional Space
									
//=========================== BSF Types =========================
struct PT_bsf_data_T {	// Data for workers
	PP_FLOAT_POINT_TYPE approximation[PP_N];	// Current approximation
};

struct PT_bsf_mapElem_T {			// Element of map list
	int rowNo;						// Row number in matrix Alpha
};

struct PT_bsf_reduceElem_T {			// Element of reduce list	
	PP_FLOAT_POINT_TYPE coordinate;		// Coordinate
};