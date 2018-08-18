/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Jacobi method implementation using BSF-M Skeleton
Module: Problem-Data.h (Problem Data)
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/

#pragma once
#include "Problem-Parameters.h"		// Problem Parameters 

//========================== Problem variables ====================================

//========================== Problem structures ====================================
static PP_FLOAT_POINT_TYPE PD_A[PP_N][PP_N];		// Coefficients of equations
static PP_FLOAT_POINT_TYPE PD_b[PP_N];			    // Vector of right parts
static PP_FLOAT_POINT_TYPE PD_Alpha[PP_N][PP_N];	// Reduced matrix
static PP_FLOAT_POINT_TYPE PD_beta[PP_N];			//  Reduced vector of right parts