#ifndef BOUND_H
#define BOUND_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_limbs.h>
#include <NTL/ZZX.h>

#include <NTL/mat_ZZ_p.h>


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

NTL::RR rhoBound(const int n, const int w, const NTL::RR& NthRootOfP);

NTL::RR nthRoot(const NTL::ZZ& P, const int n);

int calcPhi(const NTL::ZZ& Rho, const NTL::ZZX& E);

NTL::ZZ calcRho(const NTL::ZZX& M, const NTL::Mat<NTL::ZZ>& Base, const NTL::ZZX& E);

bool testPHIBound(const NTL::ZZ& P, const NTL::ZZX& E);

bool testRHOBound(const NTL::ZZ& P, const NTL::ZZX& E);

bool testPHIBoundMat(const NTL::ZZX& E, const NTL::Mat<NTL::ZZ>& Mat);

bool testRHOBoundMat(const NTL::ZZX& E, const NTL::Mat<NTL::ZZ>& Mat);

bool testPHIBoundPoly(const NTL::ZZX& E, const NTL::ZZX& M);

bool testRHOBoundPoly(const NTL::ZZX& E, const NTL::ZZX& M);

bool testDegNthRootP(const NTL::ZZ& P, int n, int deg);

bool testRootNthRootP(const NTL::ZZ& P, int n, const NTL::ZZ_p& Root);

#endif