#ifndef FINDM_H
#define FINDM_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>


#include <NTL/mat_ZZ_p.h>

#include <NTL/LLL.h>

#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h> 
#include <NTL/ZZ_pEXFactoring.h>

#include "Base.h"
#include "Norm.h"
#include "Misc.h"

/*
	function related to the search of B and M
*/

extern const char* MorBStr[2];

enum MorB : int{
	M, B
};

bool checkMorB(const NTL::ZZX& M, const NTL::ZZX& E);


bool findMorBFromRoot(const NTL::ZZ& P, const int n, const NTL::RR& W4, const NTL::ZZ_p& Root, const NTL::ZZX& E, NTL::ZZX& M_dest, NTL::Mat<NTL::ZZ>& Base_dest);

void invMorB(const NTL::ZZ& P, const NTL::ZZX& E, const NTL::ZZX& M, 
			const NTL::Mat<NTL::ZZ>& Base,  NTL::ZZX& M_inv_dest, 
			 NTL::Mat<NTL::ZZ>& Base_inv_dest);

#endif