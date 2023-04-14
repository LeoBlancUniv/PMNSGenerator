#ifndef NORM_H
#define NORM_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <NTL/mat_ZZ_p.h>

#include "Poly.h"

/*
	function related to norm calculation
*/

NTL::ZZ norm1_M_Wrapper(const NTL::ZZX& M, const NTL::ZZX& E);

NTL::ZZ norm1(const NTL::Mat<NTL::ZZ>& Mat);

NTL::ZZ norm1_M_Xn_LAM(const NTL::ZZX& M, const NTL::ZZX& E);

NTL::ZZ norm1_M_Generic(const NTL::ZZX& M, const NTL::ZZX& E);

NTL::RR calcW_Wrapper(const NTL::ZZX& E);

NTL::RR calcW_Xn_Lam(const NTL::ZZX& E);

NTL::RR calcW_Xn_X_1(const NTL::ZZX& E);

NTL::RR calcW_Xn_X2_1(const NTL::ZZX& E);

#endif