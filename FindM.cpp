#include "FindM.h"

using namespace std;
using namespace NTL;

const char* MorBStr[2] = {"M", "B"}; 

bool resultantCheck(const ZZX& M, const ZZX& E){
	/*
		compute the resultant of E and M and check that it's odd
	*/
	ZZ Val;
	ZZX U, V;

	XGCD(Val, U, V, M, E);
	bool ok = Val%2 == 1;
	
	return ok;
}

bool gcd_M_Xn_1(const ZZX& M, int n){
	/*
		compute the GCD of M viewed in F2 with X^n - 1 and check thats it's 1
	*/

	ZZ_pX Xn_p1;
	SetCoeff(Xn_p1, n);
	SetCoeff(Xn_p1, 0, 1);

	{ //in this scope, operation are done mod 2
		ZZ_pPush push;
		ZZ_p::init(ZZ(2));

		ZZ_pX res;

		res = GCD(conv<ZZ_pX>(M), Xn_p1);
		return deg(res) == 0 and IsOne(ConstTerm(res));
	}

}

bool check_Xn_Lam(const ZZX& M, const ZZX& E){
	int n = deg(E);
	ZZ Lam = -ConstTerm(E);

	bool is_pair;
	is_pair = !IsOdd(Lam);

	if (is_pair){
		//Lam is pair, we need to check if the constant terme if M is odd
		return IsOdd(ConstTerm(M));
	}
	else{
		//Lam is odd we need to check that the gcd of M and Xn_p1 is 1
		return gcd_M_Xn_1(M, n);
	}
}

bool checkMorB(const ZZX& M, const ZZX& E){
	/*
		check wether a given M is invertible mod (E, PHI)
	*/

	if (is_irreducible(E)){
		switch (_Poly_Type){
			case Xn_Lam : {

				return check_Xn_Lam(M, E);

			}break;

			default :
				return resultantCheck(M, E);
			break;
		}	
	}
	else{
		switch (_Poly_Type){
			default:
				return resultantCheck(M, E);
			break;	
		}
	}

	
}


bool findMorBFromRoot(const ZZ& P, const int n, const RR& W4, const ZZ_p& Root, const ZZX& E, ZZX& M_dest, mat_ZZ& Base_dest){
	/*
		from a root of the poly, create the associated lattice
		and tries to find depeding on _B_or_M
		-a suitable Base B
		-a suitable Poly M 
		return true or false depending on if it was found
		
		
	*/


	ZZX M;

	mat_ZZ Base = createReducedBaseFpLLL(P, n, Root);

	if (_B_or_M == MorB::B and W4 * conv<RR>(norm1(Base)) * (_Delta + 1) * (_Delta + 1) < PHI_RR){
		Base_dest = Base;
		return true;
	}

	int nbRow = Base.NumRows();
	
	ZZ Lam = -ConstTerm(E);

	int i = 0;

	RR MinBoundM = conv<RR>(0);	
	while(true){

		M = conv<ZZX>(Base[i]); //convert the ith line of the matrix into a poly

		if (checkMorB(M, E)){ //is M invertible
			RR BoundTmp = W4 * conv<RR>(norm1_M_Wrapper(M, E)) * (_Delta + 1) * (_Delta + 1);
			if (BoundTmp < PHI_RR){
				
				M_dest = M;
				Base_dest = Base;
				
				return true;
			}
			else{
				if (MinBoundM == 0){
					MinBoundM = BoundTmp;
				}
				else
				{
					if (BoundTmp < MinBoundM){
						MinBoundM = BoundTmp;
					}
				}
			}
		}

		i++;

		if (i == nbRow){ //every poly of the base was not good
			//cout << Base << endl;
			/*RR MinBound = conv<RR>(0);
			if (not IsOdd(ConstTerm(E))){
				//trying to make the sum of lines of the base work :
				ZZX A, B;
				

				for (int j = 0; j < nbRow; j++){

					bool invA;

					A = conv<ZZX>(Base[j]);

					invA = checkMorB(A, E);

					RR BoundA = W4 * conv<RR>(norm1_M_Wrapper(A, E)) * (_Delta + 1) * (_Delta + 1);
					
					

					RR BoundAB = BoundA;
					RR BoundTmpAB;
					RR BoundTmpA_B;
					RR BoundTmpB_A;

					for (int k = 0;k < nbRow; k++){
						
						bool invB;

						if (j != k){

							B = conv<ZZX>(Base[k]);

							invB = checkMorB(B, E);

							ZZX AB = A + B;

							if (checkMorB(AB, E)){
								BoundTmpAB = W4 * conv<RR>(norm1_M_Wrapper(AB, E)) * (_Delta + 1) * (_Delta + 1);
								if (MinBound == 0){
									MinBound = BoundTmpAB;
								}
								else{
									if (BoundTmpAB < MinBound){
										MinBound = BoundTmpAB;
									}
								}
							}

							ZZX A_B = A - B;
							if (checkMorB(A_B, E)){
								BoundTmpA_B = W4 * conv<RR>(norm1_M_Wrapper(A_B, E)) * (_Delta + 1) * (_Delta + 1);
								if (MinBound == 0){
									MinBound = BoundTmpA_B;
								}
								else{
									if (BoundTmpA_B < MinBound){
										MinBound = BoundTmpA_B;
									}
								}
							}

							ZZX B_A = B - A;
							if (checkMorB(B_A, E)){
								BoundTmpB_A = W4 * conv<RR>(norm1_M_Wrapper(B_A, E)) * (_Delta + 1) * (_Delta + 1);
								if (MinBound == 0){
									MinBound = BoundTmpB_A;
								}
								else{
									if (BoundTmpB_A < MinBound){
										MinBound = BoundTmpB_A;
									}
								}
							}
						}
					}
				}

				RR RR_2 = conv<RR>(2);
				RR bottom = log(RR_2);

				if (MinBoundM > MinBound){
					if (log(MinBoundM)/bottom > 80){
						//cout << Base << endl;
					}

					cout << log(MinBoundM)/bottom << " " << flush;
					if (MinBound != 0){
						cout << log(MinBound)/bottom << endl;
					}
					else{
						cout << endl;
					}
				}

				
			
				
			}*/
			return false;
		}
	}
}

void invMorB(const ZZ& P, const ZZX& E, const ZZX& M, 
			const Mat<ZZ>& Base,  ZZX& M_inv_dest,  Mat<ZZ>& Base_inv_dest){

	if (_B_or_M == MorB::B){
		Mat<ZZ> Base_inv;
		ZZ d;
		inv(d, Base_inv, Base);

		{	//create a scope where the modulus is Phi
			ZZ_pPush push;
			ZZ_p::init(PHI);

			for (int j = 0; j < Base_inv.NumRows(); j++){
				for (int k = 0; k < Base_inv.NumRows(); k++){
					ZZ_p Tmp =  conv<ZZ_p>(Base_inv[j][k]);
					Base_inv[j][k] = conv<ZZ>(inv(conv<ZZ_p>(P)) * Tmp);
				}
			}
		}

		Base_inv_dest = Base_inv;
		
	}
	else{
		//we have M here so we can calc it's inverse
		ZZX M_inv;
		ZZ Val;
		ZZX U, V;

		XGCD(Val, U, V, M, E);

		{	//create a scope where the modulus is Phi
			ZZ_pPush push;
			ZZ_p::init(PHI);

			M_inv = conv <ZZX>(conv<ZZ_pX>(U) * inv(conv<ZZ_p>(Val)));
			for (int j = 0; j < deg(M_inv) +1; j++){
				
				ZZ_p Tmp = -(conv<ZZ_p>(coeff(M_inv, j)));

				if (conv<ZZ>(Tmp) >= (PHI >>1)){
					SetCoeff(M_inv, j, conv<ZZ>(conv<ZZ>(Tmp) - PHI));
				}
				else{
					SetCoeff(M_inv, j, conv<ZZ>(Tmp));
				}
			}
		}

		M_inv_dest = M_inv;
	}
}