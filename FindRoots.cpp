#include "FindRoots.h"

using namespace std;
using namespace NTL;

const char* RootStrategyStr[6] = {
	//generic
	"Generic_Factor",
	"Generic_GCD_Xp_X",
	//Xn_Lam
	"Xn_Lam_EasyRoot",
	"Xn_Lam_NthResidue",
	"Xn_Lam_NthResidueOrFactor",
	"Xn_Lam_CompositeGCD",
	//Xn_X_1
};

//generic

bool findRootsGcdXp_X(const ZZ& P, const ZZ_pX& Poly, Vec<ZZ_p>& Roots_dest){
	/*
		generic way of finding the roots of <Poly>
		builds the poly Xp_X = X^P - X mod <Poly> and then takes R = GCD of Xp_X and <Poly>
		since P is prime, every a in Fp is a solution of Xp_X(x) = 0 mod P
		so Xp_X(x) = (x - 0)(x - 1)...(x - P-1)
		
		from this, R can be equals to : 
		- R(x) = 1, which means there where no roots in <Poly>
		- R(x) = (x - gamma_0)(x - gamma_1)...(x - gamma_b) where the gamma_i are roots of <Poly>

		in the 2nd case we find the roots by calling NTL's findRoot which gives us one of those roots 
		we then divide R(x) by (x - the root we just found), then we try to find a root again 
		repeat until R(x) is 1

		return true or false depending on if at least one roots was found

		will skip everything and return false depending on the state of _FACTORLIMIT

	*/

	if (_FactorLimit == 0){
		if (deg(Poly) >= FACTOR_LIMIT_DEFAULT){
			return false;
		}
	}
	else if (_FactorLimit > 0){
		if (deg(Poly) >= _FactorLimit){
			return false;
		}
	}



	ZZ_pX X;
	ZZ_pX Xp_X;
	
	SetCoeff(X, 1, 1);

	ZZ_pXModulus Mod(Poly);

	PowerXMod(Xp_X, P, Mod);

	sub(Xp_X, Xp_X ,X);

	ZZ_pX R = GCD(Xp_X, Poly);

	Vec<ZZ_p> Roots_vec;

	while(deg(R) != 0 or ConstTerm(R) != conv<ZZ_p>(1)){
		ZZ_p Root;
		FindRoot(Root, R);

		Roots_vec.append(Root);

		ZZ_pX D;
		SetCoeff(D, 1, 1);
		SetCoeff(D, 0, -Root);
		div(R, R, D);
	}

	if (Roots_vec.length() > 0){
		Roots_dest = Roots_vec;
		return true;
	}

	return false;

}

bool findRootsFactor(const ZZX& poly, Vec<ZZ_p>& Roots_dest){
	/*
		generic way to finds all the roots of <poly> 
		puts them in <Roots_dest> if any are found
		return true or false depending on if any roots could be found

		!!! findRootsGcdXp_X should be used instead of this !!! 

		will skip everything and return false depending on the state of _FACTORLIMIT

	*/

	if (_FactorLimit == 0){
		if (deg(poly) >= FACTOR_LIMIT_DEFAULT){
			return false;
		}
	}
	else if (_FactorLimit > 0){
		if (deg(poly) >= _FactorLimit){
			return false;
		}
	}

	
	vec_pair_ZZ_pX_long Factor_vec = berlekamp(conv<ZZ_pX>(poly));

	Vec<ZZ_p> Roots_vec;

	int nbFactor = Factor_vec.length();

	int degree = 0;

	ZZ_p Root;
	for (int i = 0; i < nbFactor; i++){

		degree = deg(Factor_vec[i].a);

		if (degree == 1){ //we found a factor that looks like (X + a), so the degree is 1
			Root = -coeff(Factor_vec[i].a, 0); //the root is minus the constant term
			Roots_vec.append(Root);
		}
	}

	if (Roots_vec.length() > 0){
		Roots_dest = Roots_vec;
		return true;
	}
	return false;
}

bool filterRoot(const ZZ& P, const ZZX& Poly, const Vec<ZZ_p>& Roots_src, Vec<ZZ_p>& Roots_dest){
	/*
		we only keep the root of the factor that pass the bound check

		if _FastRootCheck is true, will stop when one root pass the test
	*/
	int n = deg(Poly);

	Vec<ZZ_p> RootsFilter_vec;

	vec_pair_ZZX_long Factor_vec;

	Vec<ZZ_pX> FactorFilter_vec;

	factorZ(Factor_vec, Poly);

	//filter which factor or Poly are going to be considered 
	for (int i = 0; i < Factor_vec.length(); i++){
		if (testDegNthRootP(P, n, deg(Factor_vec[i].a))){
			FactorFilter_vec.append(conv<ZZ_pX>(Factor_vec[i].a));
		}
	}

	//we only keep the roots of Poly that are roots of the considered factor and that are bigger than P^(1/n) 
	for (int i = 0; i < Roots_src.length(); i++){
		bool ok = false;
		for (int j = 0; j < FactorFilter_vec.length(); j++){
			if (IsZero(eval(FactorFilter_vec[j], Roots_src[i]))){
				if (testRootNthRootP(P, n, Roots_src[i]))
					ok = true;
			}
		}
		if (ok){
			RootsFilter_vec.append(Roots_src[i]);

			if (! _DisableFastRootCheck){
				Roots_dest = RootsFilter_vec;
				return true;
			}
		}
	}

	Roots_dest = RootsFilter_vec;

	return Roots_dest.length() > 0;
}

//Xn_Lam

RootStrategy findRootStrategy_Xn_Lam(const ZZ& P, const ZZ& N, ZZ& Expo_dest, bool forceFactor){
	/*
		Find the best strategy to find roots when the poly is on the form Xn - Lam
		Expo is P-1 / B^k

	*/

	ZZ Expo = P-1;

	if (GCD (Expo, N) != 1 and GCD(Expo, N) != N){
		return RootStrategy::Xn_Lam_CompositeGCD;
	}

	int count = 0;

	while(GCD(Expo, N) == N){
		Expo /= N;
	}

	Expo_dest = Expo;

	if (count > 1){
		if (forceFactor)
			return RootStrategy::Xn_Lam_NthResidueOrFactor;
		else
			return RootStrategy::Xn_Lam_NthResidue;
	}
	else{
		return RootStrategy::Xn_Lam_EasyRoot;
	}

	/*ZZ Expo;
	ZZ B = GCD(P-1, N);
	Expo =  (P-1) / B;

	if (B == 1){ 
		Expo_dest = Expo;
		return RootStrategy::Xn_Lam_EasyRoot;
	}
	if (B == N){
		//still not convinced this works properly

		if (GCD(Expo, B) != 1){
			
			while(GCD(Expo, B) != 1)
				Expo /= B;

			Expo_dest = Expo;

			if (forceFactor)
				return RootStrategy::Xn_Lam_NthResidueOrFactor;
			else
				return RootStrategy::Xn_Lam_NthResidue;

		}
		else{

			
			
			Expo_dest = Expo;
			return RootStrategy::Xn_Lam_EasyRoot;
		}
	}
	else{
		return RootStrategy::Xn_Lam_CompositeGCD;
	}*/

}

bool findRoots_Xn_Lam(const ZZ_pX& Poly, const RootStrategy RStrat, 
	const ZZ& P, const ZZ& N, const ZZ& PoverBK, Vec<ZZ_p>& Roots_dest){

	/*
		finds the roots of <Poly> of the form Xn - Lam
		uses a lot of tricks to avoid high degree factoring
	*/


	ZZ_p Lam = -ConstTerm(Poly);
	ZZ B = GCD(P-1, N);

	switch (RStrat){
		case Xn_Lam_CompositeGCD : {
			Vec<ZZ_p> RootInit_vec;

			RootInit_vec.append(Lam);

			return findRootsGcdComposite_Xn_Lam(P, N, B, RootInit_vec, Roots_dest);
			break;
		}

		case Xn_Lam_EasyRoot : {
			//here it's guarantied that findEasyRoot works
			ZZ_p R1;
			findEasyRoot_Xn_lam(RStrat, Lam, N, PoverBK, R1);
			generateRoot_Xn_lam(P, B, R1, Roots_dest);

			return true;
			break;
		}

		case Xn_Lam_NthResidue : {

			//here it can happen that findEasyRoot doesn't return a root
			ZZ_p R1;
			if (findEasyRoot_Xn_lam(RStrat, Lam, N, PoverBK, R1)){
				generateRoot_Xn_lam(P, B, R1, Roots_dest);
				return true;
			}
			else{
				return false;
			}
			break;
		}

		case Xn_Lam_NthResidueOrFactor : {
			ZZ_p R1;
			if (findEasyRoot_Xn_lam(RStrat, Lam, N, PoverBK, R1)){
				generateRoot_Xn_lam(P, B, R1, Roots_dest);
				return true;
			}
			else{
				return findRootsGcdXp_X(P, Poly, Roots_dest);
			}
			break;
		}

		case Generic_GCD_Xp_X :
			return findRootsGcdXp_X(P, Poly, Roots_dest);
		break;

		default : cout << "ERROR in RStrat" << endl; break;

	}

	return false;
}

bool findEasyRoot_Xn_lam(const RootStrategy& RStrat, const ZZ_p& Lam, const  ZZ& N, const ZZ& Expo, ZZ_p& Roots_dest){
	/*
		tries to find a root by looking at the easy root formula
		if Lam^(Expo) = 1 mod P then Lam^u is a root where u is the inverse of N mod (Expo)
		called if RStrat is Xn_Lam_EasyRoot or Xn_Lam_NthResidue
	*/
	ZZ_p N_p = conv<ZZ_p>(N); 


	//these two pretty much do the same thing maybe factor code here
	if (RStrat == RootStrategy::Xn_Lam_EasyRoot){
		ZZ_p U;
		{ 
			
			ZZ_pPush push;
			ZZ_p::init(Expo);
			U = inv(N_p);
		}

		Roots_dest = power(Lam, conv<ZZ>(U));

		return true;
	}
	else{

		if (is_nThResidue(Lam, Expo)){
			
			ZZ_p U;
			{ 
				
				ZZ_pPush push;
				ZZ_p::init(Expo);
				U = inv(N_p);
			}

			Roots_dest = power(Lam, conv<ZZ>(U));

			return true;
		}
		else{
			return false;
		}
	}
}

bool findRootsGcdComposite_Xn_Lam(const ZZ& P, const ZZ& N, const ZZ& B, Vec<ZZ_p>& Roots_src, Vec<ZZ_p>& Roots_dest){
	/*
		only called when <B> is composite
		<B> is gcd of P-1 and <N>

		this function for each root in <Roots_src> will find the solutions of X^<N> - root
		by considering N = B * Q(<- is going to be N / B)

		looking for Nth root of A is the same as looking at the 
		Qth root of the Bth root of A

		since both B and Q are smaller than N we can hope that it makes calculation easier
		even in the worst case where we will have to use factor
		to find the root of the subpart we call findRoots_Xn_Lam

		which itself which may end up calling this function again
		on the first iteration of this function <Roots_src> will just contains the initial lambda
	*/

	//generic
	int nBRoots = Roots_src.length();


	//left part

	ZZ Left_N = B;
	//P stay the same on left and right

	Vec<ZZ_p> Left_Roots_vec; //this will hold all the Bth root

	for (int i = 0; i < nBRoots; i++){

		//for each root we try to find their Bth root

		ZZ Left_PoverBK;

		RootStrategy Left_RStrat = findRootStrategy_Xn_Lam(P, Left_N, Left_PoverBK, true);

		ZZ_pX Left_Poly;

		SetCoeff(Left_Poly, conv<int>(Left_N));
		SetCoeff(Left_Poly, 0, -Roots_src[i]);

		Vec<ZZ_p> Left_Roots_Tmp_vec;
		if (findRoots_Xn_Lam(Left_Poly, Left_RStrat, P, Left_N, Left_PoverBK, Left_Roots_Tmp_vec)){
			
			int nBRoots_tmp = Left_Roots_Tmp_vec.length();
			for (int j = 0; j < nBRoots_tmp; j++){
				Left_Roots_vec.append(Left_Roots_Tmp_vec[j]);
			}
		}

	}

	//right part, we take the Bth root and try to find their Qth root

	ZZ Right_N = N / B;
	//P stay the same on left and right
	ZZ Right_B = GCD(P-1, Right_N);

	Vec<ZZ_p> Right_Roots_vec; //this will hold all the Qth root of the Bth root

	int Left_nBRoots = Left_Roots_vec.length(); //we will iter of the Bth root

	for (int i = 0; i < Left_nBRoots; i++){
		ZZ Right_PoverBK;

		RootStrategy Right_RStrat = findRootStrategy_Xn_Lam(P, Right_N, Right_PoverBK, true);

		ZZ_pX Right_Poly;

		SetCoeff(Right_Poly, conv<int>(Right_N));
		SetCoeff(Right_Poly, 0, -Left_Roots_vec[i]);

		Vec<ZZ_p> Right_Roots_Tmp_vec; //return vec, contains all the Nth (Qth of Bth) root
		if (findRoots_Xn_Lam(Right_Poly, Right_RStrat, P, Right_N, Right_PoverBK, Right_Roots_Tmp_vec)){
			
			int nBRoots_tmp = Right_Roots_Tmp_vec.length();
			for (int j = 0; j < nBRoots_tmp; j++){
				Right_Roots_vec.append(Right_Roots_Tmp_vec[j]);
			}
		}

	}

	if (Right_Roots_vec.length() > 0){
		Roots_dest = Right_Roots_vec;
		return true;
	}
	else{
		return false;
	}

}

void generateRoot_Xn_lam(const ZZ& P, const ZZ& B, const ZZ_p& Gamma, Vec<ZZ_p>& Roots_dest){
	/*
		for XN-lam only
		build the full list of Roots from a single one
		if B = 1, there is no more roots to find (and it'll cause an infinite loop if we try)

	*/

	Vec<ZZ_p> Roots_vec;

	if (B == 1){
		Roots_vec.append(Gamma);
		Roots_dest = Roots_vec;
		return;
	}

	ZZ_p NonResidue = generateNonNthResidue(P, B);

	ZZ_p UnityRoot = power(NonResidue, (P-1)/B);

	

	for (int i = 1; i < B + 1; i++){
		Roots_vec.append(Gamma * power(UnityRoot, i));
	}

	Roots_dest = Roots_vec;

}



ZZ_p generateNonNthResidue(const ZZ& P, const ZZ& B){
	/*
		generate a non nth residue ie
		a number A so that A^(P-1/B) != 1 mod P
	*/
	ZZ_p Ret;

	ZZ Expo = (P-1) / B;

	Ret = 1; //always a nth residue

	while(is_nThResidue(Ret, Expo)){
		Ret = conv<ZZ_p>(RandomBnd(P));
	}

	return Ret;
}

bool is_nThResidue(const ZZ_p& X, const ZZ& Expo){
	return power(X, Expo) == 1;
}

RootStrategy findRootStrategy_Xn_X_1(){
	return RootStrategy::Generic_GCD_Xp_X;
}

bool findRoots_Xn_X_1(const ZZ& P, const ZZX& Poly, const RootStrategy RStrat, Vec<ZZ_p>& Roots_dest){

	switch (RStrat){
		case Generic_Factor : return findRootsFactor(Poly, Roots_dest); break;
		case Generic_GCD_Xp_X : return findRootsGcdXp_X(P, conv<ZZ_pX>(Poly), Roots_dest); break;

		default : cerr << "ERROR in RStrat" << endl; break;
	}

	return false;
}

RootStrategy findRootStrategy_Xn_X2_1(){
	return RootStrategy::Generic_GCD_Xp_X;
}

bool findRoots_Xn_X2_1(const ZZ& P, const ZZX& Poly, const RootStrategy RStrat, Vec<ZZ_p>& Roots_dest){

	switch (RStrat){
		case Generic_Factor : return findRootsFactor(Poly, Roots_dest); break;
		case Generic_GCD_Xp_X : return findRootsGcdXp_X(P, conv<ZZ_pX>(Poly), Roots_dest); break;

		default : cerr << "ERROR in RStrat" << endl; break;
	}

	return false;
}