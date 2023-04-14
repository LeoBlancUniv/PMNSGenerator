#ifndef BASE_H
#define BASE_H

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


#include <fplll/fplll.h>

#include "Misc.h"
#include "Bound.h"

/*
	Base related function
*/


NTL::Mat<NTL::ZZ> createBase(const NTL::ZZ& P, const int n, const NTL::ZZ_p& Root);

//deprecated, uses NTL's LLL alg
NTL::Mat<NTL::ZZ> createReducedBase(const NTL::ZZ& P, const int n, const NTL::ZZ_p& Root);


NTL::Mat<NTL::ZZ> createReducedBaseFpLLL(const NTL::ZZ& P, const int n, const NTL::ZZ_p& Root);

#endif
