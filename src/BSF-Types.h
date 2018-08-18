/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-M Skeleton
Module: BSF-Types.h (BSF Types)
Prefix: BT
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/
#pragma once
#include "BSF-Include.h"

struct BT_order_T {  
	char exit;		// true, if worker must stop
	int index;		// index of the first MAP sublist element
	int length;		// length of MAP sublist
	PT_bsf_data_T data;
};