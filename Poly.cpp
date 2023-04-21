#include "Poly.h"

using namespace std;
using namespace NTL;

const char* PolyTypeStr[4] = {"Xn_Lam", "Xn_X_1", "Xn_X2_1", "ANY"};

bool is_irreducible(const ZZX& E){
	/*
		tries to factor the poly <E> to see if it is irreducible
		generic
	*/
	vec_pair_ZZX_long Factor_vec;
	ZZ C;
	factor(C, Factor_vec, E);
	//cout << Factor_vec << endl;
	return Factor_vec.length() == 1;
}

void factorZ(vec_pair_ZZX_long& Factor_vec_dest, const ZZX& Poly){
	/*
		returns the list of factor of <Poly>
		generic
	*/
	vec_pair_ZZX_long Factor_vec;
	ZZ C;
	factor(C, Factor_vec, Poly);
	Factor_vec_dest = Factor_vec;
}

bool polyInitWrapper(int n, ZZX& Poly){

	/*
		wrapper to initialize <Poly> with the form _Poly_Type
		return false if the poly can't be initialized
	*/

	switch (_Poly_Type){
		case Xn_Lam : return polyInit_Xn_Lam(n, Poly); break;
		case Xn_X_1 : return polyInit_Xn_X_1(n, Poly); break;
		case Xn_X2_1 : return polyInit_Xn_X2_1(n, Poly); break;
		default : cout << "ERROR in _Poly_Type" << endl; break;
	}

	return false;

}

bool polyInit_Xn_Lam(int n, ZZX& Poly){
	/*
		<Poly> init for Xn_Lam
	*/

	ZZX Ret;
	int lam = 0; //to remove unitialized warning

	if (_Xn_Lam_Lam_Set){
		lam = _Xn_Lam_Lam;
	}

	if (_Xn_Lam_Lam_RangeMin != -1){
		lam = _Xn_Lam_Lam_RangeMin;
	}else{
		lam = LAM_START_DEFAULT;
	}

	SetCoeff(Ret, n);
	SetCoeff(Ret, 0, -lam);

	Poly = Ret;
	return true;
}

bool polyInit_Xn_X_1(int n, ZZX& Poly){
	/*
		<Poly> init for Xn_X_1
	*/
	ZZX Ret;

	SetCoeff(Ret, n);
	SetCoeff(Ret, 1, 1);
	SetCoeff(Ret, 0, 1);

	Poly = Ret;

	return true;
}

bool polyInit_Xn_X2_1(int n, ZZX& Poly){
	/*
		<Poly> init for Xn_X2_1
	*/
	if (n%2 == 1){
		return false;
	}

	ZZX Ret;

	SetCoeff(Ret, n);
	SetCoeff(Ret, n/2, 1);
	SetCoeff(Ret, 0, 1);

	Poly = Ret;

	return true;
}

//

bool polyNextWrapper(ZZX& Poly_dest, const ZZX& Poly_src){
	/*
		wrapper for finding from <Poly_src> the next poly to check
		return false if the next poly can't be made 
	*/
	switch (_Poly_Type){
		case Xn_Lam : return polyNext_Xn_Lam(Poly_dest, Poly_src); break;
		case Xn_X_1 : return polyNext_Xn_X_1(Poly_dest, Poly_src); break;
		case Xn_X2_1 : return polyNext_Xn_X2_1(Poly_dest, Poly_src); break;
		default : cout << "ERROR in _Poly_Type" << endl; break;
	}
	return false;
}

bool polyNext_Xn_Lam(ZZX& Poly_dest, const ZZX& Poly_src){
	/*
		Xn_Lam
		if lam is pos -> lam = -lam
		if lam is neg -> lam = -lam + 1

		will respect the min and max range on lam if set in meta.h
	*/

	int n = deg(Poly_src);
	int lam = conv<int>(-ConstTerm(Poly_src));

	if (_Xn_Lam_Lam_Set){
		return false;
	}

	if (lam > 0){
		lam = - lam;
	}
	else{
		lam = -lam;
		lam++;
		if (_Xn_Lam_Lam_RangeMax != -1){
			if (lam > _Xn_Lam_Lam_RangeMax){
				return false;
			}
		}
	}

	ZZX Ret;

	SetCoeff(Ret, n);
	SetCoeff(Ret, 0, -lam);

	Poly_dest = Ret;
	return true;


}

bool polyNext_Xn_X_1(ZZX& Poly_dest, const ZZX& Poly_src){
	/*
		Xn_X_1
		makes variation on the sign in front of the X and const term
	*/

	ZZ c, x;

	int n = deg(Poly_src);

	if (ConstTerm(Poly_src) > 0){
		c = -1;
		x = coeff(Poly_src, 1);
	}
	else if (coeff(Poly_src, 1) < 0){
		return false;
	}
	else{
		c = 1;
		x = -1;
	}


	ZZX Ret;

	SetCoeff(Ret, n);
	SetCoeff(Ret, 1, x);
	SetCoeff(Ret, 0, c);

	Poly_dest = Ret;

	return true;
}

bool polyNext_Xn_X2_1(ZZX& Poly_dest, const ZZX& Poly_src){
	/*
		Xn_X2_1
		makes variation on the X2 term
	*/
	int n = deg(Poly_src);

	ZZ N2 = coeff(Poly_src, n/2);

	if (N2 < 0){
		return false;
	}

	ZZX Ret;

	SetCoeff(Ret, n);
	SetCoeff(Ret, n/2, -1);
	SetCoeff(Ret, 0, 1);

	Poly_dest = Ret;

	return true;

}

//

bool polyCheckWrapper(const ZZX& Poly, const ZZ& P){

	/*
		wrapper for check function 
		those function should return false if the poly can't be used , true if they can
	*/

	switch (_Poly_Type){
		case Xn_Lam : return polyCheck_Xn_Lam(Poly, P); break;
		case Xn_X_1 : return polyCheck_Xn_X_1(Poly); break;
		case Xn_X2_1 : return polyCheck_Xn_X2_1(Poly); break;
		default : cout << "ERROR in _Poly_Type" << endl; break;
	}
	return false;
}

bool polyCheck_Xn_Lam(const ZZX& Poly, const ZZ& P){
	/*
		fast check if <Poly> contains root
	*/

	ZZ N = conv<ZZ>(deg(Poly));

	ZZ B = GCD(P-1, N);

	ZZ_p Lam = conv<ZZ_p>(-ConstTerm(Poly));

	return is_nThResidue(Lam, (P-1)/B);
}



bool polyCheck_Xn_X_1(const ZZX& Poly){
	//if you find a condidtion to use the Xn_X_1 it goes here
	return true;
}

bool polyCheck_Xn_X2_1(const ZZX& Poly){
	//if you find a condidtion to use the Xn_X2_1 it goes here

	return true;
}

bool polyNothingFoundWrapper(RootStrategy& RStrat){

	/*
		wrapper for when no pmns was found for a given n
		you can update the RStrat
	*/

	switch (_Poly_Type){
		case Xn_Lam: return polyNothingFound_Xn_Lam(RStrat); break;
		case Xn_X_1: return polyNothingFound_Xn_X_1(RStrat); break;
		case Xn_X2_1: return polyNothingFound_Xn_X2_1(RStrat); break;
		default : cout << "ERROR in _Poly_Type" << endl; break;
	}

	return false;
}

bool polyNothingFound_Xn_Lam(RootStrategy& RStrat){
	/*
		Xn_Lam
		NthResidue not finding roots doesn't guaranties that there were no root to find
		so we check again in GCDXp_X mode
	*/
	if (RStrat == RootStrategy::Xn_Lam_NthResidue){
		RStrat = RootStrategy::Generic_GCD_Xp_X;
		return false;
	}
	else{
		return true;
	}
}

bool polyNothingFound_Xn_X_1(RootStrategy& RStrat){
	return true;
}

bool polyNothingFound_Xn_X2_1(RootStrategy& RStrat){
	return true;
}

void polyPrintWrapper(const ZZX& Poly){
	/*
		_Verbose print wrapper
	*/
	switch (_Poly_Type){
		case Xn_Lam : polyPrint_Xn_Lam(Poly); break;
		case Xn_X_1 : polyPrint_Xn_X_1(Poly); break;
		case Xn_X2_1 : polyPrint_Xn_X2_1(Poly); break;
		default : cout << "ERROR in _Poly_Type" << endl; break;
	}


}

void polyPrint_Xn_Lam(const ZZX& Poly){
	int n = deg(Poly);
	int lambda = conv<int>(-ConstTerm(Poly));


	cout << "iter (n = " << n << ", lambda = " << lambda << ")" << endl;
}

void polyPrint_Xn_X_1(const ZZX& Poly){
	ZZ x, c;
	int n = deg(Poly);
	x = coeff(Poly, 1);
	c = ConstTerm(Poly);

	cout << "iter (n = " << n << ", x = " << x << ", c = " << c << ")" << endl;
}

void polyPrint_Xn_X2_1(const ZZX& Poly){
	ZZ n2, c;
	int n = deg(Poly);
	n2 = coeff(Poly, n/2);
	c = ConstTerm(Poly);

	cout << "iter (n = " << n << ", n/2 = " << n2 << ", c = " << c << ")" << endl;
}