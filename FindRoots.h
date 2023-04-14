#ifndef FINDROOTS_H
#define FINDROOTS_H

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

#include <NTL/RR.h>

#include <omp.h>

#include "Misc.h"
#include "Bound.h"
#include "Meta.h"

/*
	function to find roots from poly
	contains generic method that will work on any form of poly
	and specific way that only work on particular form
*/

extern const char* RootStrategyStr[6];

enum RootStrategy : int
{
	//generic
	Generic_Factor,
	Generic_GCD_Xp_X,
	//Xn_Lam
	Xn_Lam_EasyRoot,
	Xn_Lam_NthResidue,
	Xn_Lam_NthResidueOrFactor,
	Xn_Lam_CompositeGCD,
	//Xn_X_1

	//Xn_X2_1
};


//generic

bool findRootsGcdXp_X(const NTL::ZZ& P, const NTL::ZZ_pX& Poly, NTL::Vec<NTL::ZZ_p>& Roots_dest);

bool findRootsFactor(const NTL::ZZX& Poly, NTL::Vec<NTL::ZZ_p>& Roots_dest);

bool filterRoot(const NTL::ZZ& P, const NTL::ZZX& Poly, 
				const NTL::Vec<NTL::ZZ_p>& Roots_src, NTL::Vec<NTL::ZZ_p>& Roots_dest);

//Xn_Lam

RootStrategy findRootStrategy_Xn_Lam(const NTL::ZZ& P, const NTL::ZZ& N, 
									NTL::ZZ& Expo_dest, bool forceFactor = false);

bool findEasyRoot_Xn_lam(const RootStrategy& RStrat, const NTL::ZZ_p& Lam, 
						const NTL::ZZ& N, const NTL::ZZ& Expo, NTL::ZZ_p& Roots_dest);


bool findRoots_Xn_Lam(const NTL::ZZ_pX& Poly, const RootStrategy RStrat, 
				const NTL::ZZ& P, const NTL::ZZ& N, 
				const NTL::ZZ& PoverBK, NTL::Vec<NTL::ZZ_p>& Roots_dest);

bool findRootsGcdComposite_Xn_Lam(const NTL::ZZ& P, const NTL::ZZ& N, const NTL::ZZ& B, 
					NTL::Vec<NTL::ZZ_p>& Roots_src, NTL::Vec<NTL::ZZ_p>& Roots_dest);

void generateRoot_Xn_lam(const NTL::ZZ& P, const NTL::ZZ& B, const NTL::ZZ_p& Gamma, NTL::Vec<NTL::ZZ_p>& Roots_dest);

NTL::ZZ_p generateNonNthResidue(const NTL::ZZ& P, const NTL::ZZ& B);

bool is_nThResidue(const NTL::ZZ_p& X, const NTL::ZZ& Expo);

//Xn_X_1

RootStrategy findRootStrategy_Xn_X_1();

bool findRoots_Xn_X_1(const NTL::ZZ& P, const NTL::ZZX& Poly, const RootStrategy RStrat, NTL::Vec<NTL::ZZ_p>& Roots_dest);


//Xn_X2_1

RootStrategy findRootStrategy_Xn_X2_1();

bool findRoots_Xn_X2_1(const NTL::ZZ& P, const NTL::ZZX& Poly, const RootStrategy RStrat, NTL::Vec<NTL::ZZ_p>& Roots_dest);


#endif