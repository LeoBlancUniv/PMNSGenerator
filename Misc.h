#ifndef MISC_H
#define MISC_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_limbs.h>
#include <NTL/ZZX.h>

#include <NTL/mat_ZZ_p.h>

#include <NTL/RR.h>

#include <gmp.h>
#include <fplll/fplll.h>

#include <fstream>
#include <string.h>

#include "FindRoots.h"

/*
	utilitary function that don't have their place anywhere else
*/









NTL::Vec<int> primeFactor(int n);

void convertZZtoMPZ(mpz_t Dest, const NTL::ZZ& Src);
void convertMPZtoZZ(NTL::ZZ& Dest, mpz_t Src);

void convertMatZZtoMatMPZ(fplll::ZZ_mat<mpz_t>& Dest, const NTL::Mat<NTL::ZZ>& Src);
void convertMatMPZtoMatZZ(NTL::Mat<NTL::ZZ>& Dest, const fplll::ZZ_mat<mpz_t>& Src);

void freeMatMPZ(fplll::ZZ_mat<mpz_t>& Mat);
void writeToFile(std::string filename, const NTL::ZZ& P, const int n, const NTL::ZZ_p& Gamma,
				 const NTL::ZZ& Rho, const NTL::ZZX& E, const NTL::Mat<NTL::ZZ>& Base,
				 const NTL::Mat<NTL::ZZ>& Base_inv, const NTL::ZZX& M, const NTL::ZZX& M_inv);

void writeVec(std::ofstream& File, NTL::Vec<NTL::ZZ> Vec);
void writeMatrix(std::ofstream& File, NTL::Mat<NTL::ZZ> Mat);
#endif