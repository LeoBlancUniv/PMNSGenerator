#include "Bound.h"

using namespace std;
using namespace NTL;

RR BoundCheckHack::Bound = RR(0);

BoundCheckHack::BoundCheckHack(){

}

BoundCheckHack::~BoundCheckHack(){
}

void BoundCheckHack::updateBound(const ZZ& P, const int n){
	//update the bound
	BoundCheckHack::Bound = nthRoot(P, n);
}

long BoundCheckHack::getBoundFunc(const Vec<ZZ>& Vec){
	/*
		forces the LLL to stop if we can find a vector where every coef is below the bound
	*/
	for (int i = 0; i < Vec.length(); i++){
		if (conv<RR>(abs(Vec[i])) > BoundCheckHack::Bound){
			return 0;
		}
	}
	//cout << " LLL ended early " << flush;
	return -1;
}

RR phiBound(const int n, const int w, const RR& nthRootoFp){
	/*
		this is the expected bound after LLL
	*/

	RR Tmp;

	Tmp = pow(1.02, n) * 4 * w * nthRootoFp * (_Delta + 1) * (_Delta + 1);

	if (_AugmentedMinN){
		Tmp *= nthRoot(conv<ZZ>(n), 2);
	}



	return Tmp;
}

RR nthRoot(const ZZ& P, const int n){
	/*
		returns p^1/n
	*/

	RR P_RR = conv<RR>(P);

	RR RR_2 = conv<RR>(2);

	RR up = log(P_RR);
	RR bottom = log(RR_2);

	RR Log2P = up / bottom;

	return pow(RR_2, Log2P / n);


}

ZZ calcRho(const ZZX& M, const Mat<ZZ>& Base, const ZZX& E){
	/*
		calculate Rho depending if we where looking for M or B
	*/

	ZZ Rho = ZZ(0);

	if (_B_or_M == MorB::B){
		while (ZZ(1) << conv<int>(Rho) < 2 * norm1(Base)){
			Rho++;
		}
	}
	if (_B_or_M == MorB::M){
		while (ZZ(1) << conv<int>(Rho) < 2 * norm1_M_Wrapper(M, E)){
			Rho++;
		}
	}

	return Rho;
	
}

bool testPHIBound(const ZZ& P, const ZZX& E){
	/*
		check if the poly <E> is in the PHI bound 
	*/
	int n = deg(E);
	RR W4 = 4 * calcW_Wrapper(E);
	RR NthRootOfP = nthRoot(P, n);
	RR BoundLimit = NthRootOfP * W4 * (_Delta + 1) * (_Delta + 1);

	if (_ExactMaxBound){

		//we mult the bound by the deg of the highest factor of E

		vec_pair_ZZX_long Factor_vec;

		factorZ(Factor_vec, E);

		int m = 0;

		for (int i = 0; i < Factor_vec.length(); i++){
			if (deg(Factor_vec[i].a) > m)
				m = deg(Factor_vec[i].a);
		}

		BoundLimit *= m;
	}

	return BoundLimit < PHI_RR;
}

bool testDegNthRootP(const ZZ& P, int n, int deg){
	/*
		check if deg is less than the nth root of P
	*/

	RR NthRootOfP = nthRoot(P, n);
	ZZ NthRootOfP_ZZ = conv<ZZ>(NthRootOfP);
	

	return conv<ZZ>(deg) < NthRootOfP_ZZ;
}

bool testRootNthRootP(const ZZ& P, int n, const ZZ_p& Root){

	/*
		check if the root is bigger than nth root of P
	*/

	RR NthRootOfP = nthRoot(P, n);
	ZZ NthRootOfP_ZZ = conv<ZZ>(NthRootOfP);
	ZZ Root_ZZ = abs(conv<ZZ>(Root));

	return Root_ZZ >  NthRootOfP_ZZ and Root_ZZ < P -NthRootOfP_ZZ;
}