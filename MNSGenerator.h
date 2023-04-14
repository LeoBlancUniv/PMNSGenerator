#ifndef MNSGENERATOR_H
#define MNSGENERATOR_H

#include <NTL/RR.h>

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

#include <ctime>
#include <string>

#include <omp.h>

#include "Norm.h"
#include "FindM.h"
#include "FindRoots.h"
#include "Misc.h"
#include "Meta.h"


//this generation will be optimised for E = X^n - lamda
bool testFullGenerate(const NTL::ZZ& P, int n, const NTL::ZZX& E, 
						NTL::ZZ& Rho_dest, 
						NTL::ZZ_p& Gamma_dest, NTL::mat_ZZ& Base_dest, NTL::mat_ZZ& Base_inv_dest, 
						NTL::ZZX& M_dest, NTL::ZZX& M_inv_dest);


void MnsGenerator();

void test(int order);

void MNSGenerator();


#endif