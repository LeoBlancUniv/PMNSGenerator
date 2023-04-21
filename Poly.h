#ifndef POLY_H
#define POLY_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <NTL/mat_ZZ_p.h>

#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h> 

#include "Meta.h"
#include "FindRoots.h"

/*
	function to manipulate Poly E with different form
*/

extern const char* PolyTypeStr[4]; 

enum RootStrategy : int;

enum PolyType : int
{
	Xn_Lam,
	Xn_X_1,
	Xn_X2_1,
	ANY
};

bool is_irreducible(const NTL::ZZX& E);
void factorZ(NTL::vec_pair_ZZX_long& Factor_vec_dest, const NTL::ZZX& Poly);

bool polyInitWrapper(int n, NTL::ZZX& Poly);
bool polyNextWrapper(NTL::ZZX& Poly_dest, const NTL::ZZX& Poly_src);
bool polyCheckWrapper(const NTL::ZZX& Poly, const NTL::ZZ& P);
void polyPrintWrapper(const NTL::ZZX& Poly);
bool polyNothingFoundWrapper(RootStrategy& RStrat);

bool polyInit_Xn_Lam(int n, NTL::ZZX& Poly);
bool polyNext_Xn_Lam(NTL::ZZX& Poly_dest, const NTL::ZZX& Poly_src);
bool polyCheck_Xn_Lam(const NTL::ZZX& Poly, const NTL::ZZ& P);
void polyPrint_Xn_Lam(const NTL::ZZX& Poly);
bool polyNothingFound_Xn_Lam(RootStrategy& RStrat);

bool polyInit_Xn_X_1(int n, NTL::ZZX& Poly);
bool polyNext_Xn_X_1(NTL::ZZX& Poly_dest, const NTL::ZZX& Poly_src);
bool polyCheck_Xn_X_1(const NTL::ZZX& Poly);
void polyPrint_Xn_X_1(const NTL::ZZX& Poly);
bool polyNothingFound_Xn_X_1(RootStrategy& RStrat);

bool polyInit_Xn_X2_1(int n, NTL::ZZX& Poly);
bool polyNext_Xn_X2_1(NTL::ZZX& Poly_dest, const NTL::ZZX& Poly_src);
bool polyCheck_Xn_X2_1(const NTL::ZZX& Poly);
void polyPrint_Xn_X2_1(const NTL::ZZX& Poly);
bool polyNothingFound_Xn_X2_1(RootStrategy& RStrat);








#endif