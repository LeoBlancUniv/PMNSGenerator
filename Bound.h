#ifndef BOUND_H
#define BOUND_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_limbs.h>
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

#include "Meta.h"
#include "Misc.h"
#include "Norm.h"

/*
	function related to the phi bound
*/

class BoundCheckHack{
	/*
		class to bypass a restriction for NTL's LLL check function
	*/
public:

	BoundCheckHack();
	~BoundCheckHack();

	static NTL::RR Bound;

	static void updateBound(const NTL::ZZ& P, const int n);
	static long getBoundFunc(const NTL::Vec<NTL::ZZ>& Vec);
};

NTL::RR phiBound(const int n, const int w, const NTL::RR& nthRootToFp);

NTL::RR nthRoot(const NTL::ZZ& P, const int n);

NTL::ZZ calcRho(const NTL::ZZX& M, const NTL::Mat<NTL::ZZ>& Base, const NTL::ZZX& E);

bool testPHIBound(const NTL::ZZ& P, const NTL::ZZX& E);

bool testDegNthRootP(const NTL::ZZ& P, int n, int deg);

bool testRootNthRootP(const NTL::ZZ& P, int n, const NTL::ZZ_p& Root);

#endif