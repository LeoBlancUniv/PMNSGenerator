#include "Norm.h"

using namespace std;
using namespace NTL;

NTL::ZZ norm1_M_Wrapper(const NTL::ZZX& M, const NTL::ZZX& E){

	/*
		wrapper that decide how to calculate the norm 1 of <M> depending on the form of <E>
	*/

	switch(_Poly_Type){

		case Xn_Lam :
			return norm1_M_Xn_LAM(M, E);
		break;

		default :
			return norm1_M_Generic(M, E);
		break;
	}

}

ZZ norm1(const mat_ZZ& Mat){
	/*
		calc the norm1 of <Mat>
		generic
	*/
	ZZ Sum, Tmp;

	Sum = 0;

	int nbRow = Mat.NumRows();

	int nbCol = Mat.NumCols();

	for (int i = 0; i < nbRow; i++){
		Tmp = 0;
		for (int j = 0; j < nbCol; j++){
			Tmp += abs(Mat[j][i]);
		}

		if (Tmp > Sum)
			Sum = Tmp;
	}

	return Sum;
}

ZZ norm1_M_Xn_LAM(const ZZX& M, const ZZX& E){
	/*
		optimised version of the norm1 of poly M
		only works if the poly is on the form Xn-lam
		m0 + Lam(m1 + .. + mi)

	*/

	ZZ Lam = -ConstTerm(E);

	ZZ Sum, SumTmp;

	Sum = abs(ConstTerm(M)); 

	int nbCoeff = deg(M) + 1;


	SumTmp = 0;
	for (int i = 1; i < nbCoeff; i++) { // start at one since the first coeff isn't mult by Lam
		SumTmp += abs(coeff(M, i));
	}

	SumTmp *= abs(Lam);

	return Sum + SumTmp;
}

ZZ norm1_M_Generic(const ZZX& M, const ZZX& E){
	/*
		create the matrix assosiated with M to take it's norm1
		generic
	*/


	ZZ_pX E_pX = conv<ZZ_pX>(E);

	int n = deg(E_pX);

	Mat<ZZ> M_Mat;
	M_Mat.SetDims(n, n);

	ZZX Tmp(M);

	ZZX X;
	SetCoeff(X, 1, 1);

	for (int i = 0; i < n; i++){
		Vec<ZZ> Tmp_vec = conv<Vec<ZZ>>(Tmp);
		if (Tmp_vec.length() < n){

			int padc = n - Tmp_vec.length();

			for (int j = 0; j < padc; j++){
				Tmp_vec.append(ZZ(0));
			}

			
		}
		M_Mat[i] = Tmp_vec;
		MulMod(Tmp, Tmp, X, E);


	}

	ZZ Ret = norm1((M_Mat));
	return Ret;

}

NTL::RR calcW_Wrapper(const NTL::ZZX& E){

	/*
		wrapper to decide how to calculate W depending on the form of <E>
	*/

	switch (_Poly_Type){
		case Xn_Lam : return calcW_Xn_Lam(E); break;
		case Xn_X_1 : return calcW_Xn_X_1(E); break;
		case Xn_X2_1 : return calcW_Xn_X2_1(E); break;
		default: cout << "_Poly_Type not set" << endl; exit(-1);
	}
}

NTL::RR calcW_Xn_Lam(const NTL::ZZX& E){
	/*
		W calc for Xn-lam
	*/
	int n = deg(E);
	int lam = conv<int>(-ConstTerm(E)); //to check

	return conv<RR>((abs(lam) * (n-1)) + 1);

}

NTL::RR calcW_Xn_X_1(const NTL::ZZX& E){
	/*
		W calc for Xn_X_1
	*/
	int n = deg(E);

	return conv<RR>((2 * n) - 1);
}

NTL::RR calcW_Xn_X2_1(const NTL::ZZX& E){
	/*
		W calc for Xn_X2_1
	*/
	int n = deg(E);

	return conv<RR>((3 * n) / 2);
}