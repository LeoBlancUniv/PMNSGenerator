#include "Base.h"

using namespace std;
using namespace NTL; 	
using namespace fplll;


Mat<ZZ> createBase(const ZZ& P, const int n, const ZZ_p& Root){
	/*
		create the basis of a lattice with gamma
	*/

	Mat<ZZ> Base;
	Base.SetDims(n, n);

	//we are going to fill the Base in 3 step

	//1 : Base[0][0] set to p, the rest of the row to 0 (this is done by default)
	Base[0][0] = P;

	//2 : we fill the first column (except the first entry since it's P) with -gamma^i mod n
	ZZ_p Tmp(Root);

	for (int i = 1; i < n; i++){
		Base[i][0] = -conv<ZZ>(Tmp);
		Tmp *= Root;
	}

	//3 : add identity matrix on the bottom right on the Base

	for (int i = 1; i < n; i++){
		for (int j = 1; j < n; j++){
			if (i == j)
				Base[i][j] = ZZ(1);
		}
	}
	//cout << Base << endl;
	return Base;
}

Mat<ZZ> createReducedBase(const ZZ& P, const int n, const ZZ_p& Root){
	/*
		use the LLL alg on the base
		not really usefull since the FpLLL function is faster
	*/

	

	Mat<ZZ> Base = createBase(P, n, Root);

	BoundCheckHack::updateBound(P, n);

	/*if (n < 10){
		G_LLL_FP(Base, 0.99, 0, BoundCheckHack::getBoundFunc);
	}
	else if (n < 37){
		G_LLL_FP(Base, 0.99, 0, BoundCheckHack::getBoundFunc);
	}
	else if (n < 162){
		G_LLL_XD(Base, 0.99, 0, BoundCheckHack::getBoundFunc);
	}
	else{
		G_LLL_RR(Base, 0.99, 0, BoundCheckHack::getBoundFunc);
	}*/

	
	G_BKZ_RR(Base);
	
	
	return Base;
}

Mat<ZZ> createReducedBaseFpLLL(const ZZ& P, const int n, const ZZ_p& Root){
	/*
		create reduced base from the root using FpLLL's LLL alg
		FpLLL seems to free on it's own the GMP object passed to him
	*/

	Mat<ZZ> Base = createBase(P, n, Root);
	Mat<ZZ> BaseLLL;

	ZZ_mat<mpz_t> Base_MPZ;

	//switch to something fpLLL understand
	convertMatZZtoMatMPZ(Base_MPZ, Base);



	ZZ_mat<mpz_t> U;
	ZZ_mat<mpz_t> U_inv;


	//FpLLL call

	if (_BKZ){
		Wrapper W = Wrapper(Base_MPZ, U, U_inv, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);

		W.lll();
		bkz_reduction(Base_MPZ, 10, BKZ_DEFAULT, FT_MPFR, 500);
		//bkz_reduction(Base_MPZ, 10);

	}
	else{
		Wrapper W = Wrapper(Base_MPZ, U, U_inv, LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);

		W.lll();
	}
	


	convertMatMPZtoMatZZ(BaseLLL, Base_MPZ);

	freeMatMPZ(Base_MPZ);

	return BaseLLL;
}
