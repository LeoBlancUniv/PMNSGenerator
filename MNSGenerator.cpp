#include "MNSGenerator.h"


using namespace std;
using namespace NTL;





bool testFullGenerate(const ZZ& P, int n, const ZZX& E, 
						Vec<ZZ_p>& Roots_vec, 
						ZZ& Rho_dest, ZZ_p& Gamma_dest, 
						mat_ZZ& Base_dest, mat_ZZ& Base_inv_dest,
						ZZX& M_dest, ZZX& M_inv_dest){
	/*
		from a poly E, p, and n, tries to generate a PMNS
		return true or false depending on the succes
	*/
	
	

	ZZ Rho;
	mat_ZZ Base;


	int nBRoots;

	
	nBRoots = Roots_vec.length();
	
	while(true){		

		for (int i = 0; i < nBRoots; i++){

			ZZX M;



			if (findMorBFromRoot(P, n, Roots_vec[i], E, M, Base)){
				invMorB(P, E, M, Base, M_inv_dest, Base_inv_dest);
				if (_FixedRho){
					Rho = calcRho(M, Base, E);
					_Phi_Log = calcPhi(Rho, E);

					PHI = ZZ(1) << _Phi_Log;
					PHI_RR = conv<RR>(PHI);
				}
				else{
					Rho = calcRho(M, Base, E);
				}
				
				
				Base_dest = Base;
				M_dest = M;
				Rho_dest = Rho;
				Gamma_dest = Roots_vec[i];
				return true;
			}
		}
		return false;
	}
}



bool generateFromE(const ZZ& P, const int n, const RootStrategy RStrat, const ZZ& PoverBK, const ZZX& E,
					ZZ& Rho_dest, ZZ_p& Gamma_dest, 
					mat_ZZ& Base_dest, mat_ZZ& Base_inv_dest, 
					ZZX& M_dest, ZZX& M_inv_dest){

	/*
		extract the roots of <E> and send them to the next function
	*/

	

	Vec<ZZ_p> Roots_vec;

	bool rootsFound = false;

	if (_SVerbose){
		polyPrintWrapper(E);
	}

	switch (_Poly_Type){
		case Xn_Lam :
			
			if (findRoots_Xn_Lam(conv<ZZ_pX>(E), RStrat, P, conv<ZZ>(n), PoverBK, Roots_vec)){
				rootsFound = true;
			}
			

		break;

		case Xn_X_1 :

			if (findRoots_Xn_X_1(P, E, RStrat, Roots_vec)){
				rootsFound = true;
			}

		break;

		case Xn_X2_1 :
			if (findRoots_Xn_X2_1(P, E, RStrat, Roots_vec)){
				rootsFound = true;
			}
		break;

		default : cout << "ERROR in _Poly_Type" << endl; break;
	}
	
	if (rootsFound){ 
			
		Vec<ZZ_p> RootsFilter_vec;

		if (filterRoot(P, E, Roots_vec, RootsFilter_vec)){
			//we can try to generate a pmns
			if (testFullGenerate(P, n, E, RootsFilter_vec,
										Rho_dest, Gamma_dest, Base_dest, Base_inv_dest
										, M_dest, M_inv_dest)){
				return true;
			}	
		}
	}
	
	return false;
}











bool generateFromN(ZZ& P, int n){
	RR nthRootOfP = nthRoot(P, n);

	mat_ZZ Base, Base_inv;
	ZZ_p Gamma;
	ZZ Rho;
	ZZX M, M_inv;

	ZZX TmpE;
	Vec<ZZX> E_vec;

	RootStrategy RStrat;

	//Xn_Lam
	ZZ PoverBK;

	//Xn_X_1





	switch (_Poly_Type){
		case Xn_Lam : RStrat = findRootStrategy_Xn_Lam(P, conv<ZZ>(n), PoverBK); break;
		case Xn_X_1 : RStrat = findRootStrategy_Xn_X_1(); break;
		case Xn_X2_1 : RStrat = findRootStrategy_Xn_X2_1(); break;
		default: cerr << "_Poly_Type not set" << endl; exit(-1);
	}

	if (_Verbose){
		cout << "(PolyType : " << PolyTypeStr[_Poly_Type]
			 << ", n : " << n << ", Strategy : " << RootStrategyStr[RStrat] << ") : " << flush;
	}

	//E_vec init

	if (!polyInitWrapper(n, TmpE)){
		if (_Verbose){
			cout << "No Good Poly found" << endl;
		}
		
		return false;
	}

	if (_FixedRho ? testRHOBound(P, TmpE) : testPHIBound(P, TmpE)){
		if (polyCheckWrapper(TmpE, P)){
			E_vec.append(TmpE);
		}
		
	}
	else{
		//no valid poly found, we can stop here
		if (_Verbose){
			cout << "No Good Poly found" << endl;
		}

		return false;
	}

	while((polyNextWrapper(TmpE, TmpE)) and (_FixedRho ? testRHOBound(P, TmpE) : testPHIBound(P, TmpE))){
		
		if (polyCheckWrapper(TmpE, P)){
			E_vec.append(TmpE);
		}
	}

	if (E_vec.length() == 0){
		if (_Verbose){
			cout << "No Good Poly found" << endl;
		}
		return false;
	}	

	if (_Verbose){
		cout << "Searching through " << E_vec.length() << " possible polynomials" << endl;
	}

	
	bool found = false;
	int nb;
	

	omp_set_num_threads(_Nb_Thread);
	#pragma omp parallel for
	for (int i = 0; i < E_vec.length(); i++){

		bool skip;

		ZZ_p::init(P); 

		#pragma omp flush
		#pragma omp atomic read
		skip = found; 


		if (!skip){
			if (generateFromE(P, n, RStrat, PoverBK, E_vec[i], Rho, Gamma, Base, Base_inv, M, M_inv)){

				#pragma omp flush(found)
				#pragma omp atomic write
					found = true;
				
				#pragma omp flush(nb)
				#pragma omp atomic write
					nb = i;
			}
		}
	}

	

	if (found){
		if (_Show){
			cout << "PMNS FOUND : " << endl;
			cout << "p : " << P << endl;
			cout << "n : " << n << endl;
			cout << "gamma : " << Gamma << endl;
			cout << "rho : " << Rho << endl;
			cout << "phi : " << _Phi_Log << " " << PHI << endl;
			cout << "E : " << E_vec[nb] << endl;

			if (_B_or_M == MorB::B){
				cout << "Base : " << Base << endl;
				cout << "Base_inv : " << Base_inv << endl;
			}
			else{
				cout << "M : " << M << endl;
				cout << "M_inv : " << M_inv << endl;
			}
			
		}

		if (_Out){
			writeToFile(_Path + _Filename, P, n, Gamma, Rho, E_vec[nb], Base, Base_inv, M, M_inv);			
		}
	}


	return found;
}

void MNSGenerator(){

	/*
		this function's goal is to setup all the parameter to call the generator
	*/



	if (_Verbose){
		cout << "Starting PMNS generation" << endl;
		cout << "Looking for " << MorBStr[_B_or_M] << endl;
		cout << "PHI : " << _Phi_Log << endl;
		cout << "Delta : " << _Delta << endl;
		
	}

	ZZ P;
	int nbBitsP = 0; //supprese uninitialised var warning

	if (_P != -1){
		P = _P;
		nbBitsP = NumBits(P);
	}
	if (_P_ORDER != -1){
		nbBitsP = 1 << _P_ORDER;
	}
	if (_P_Size != -1){
		nbBitsP = _P_Size;
	}
	
	if (_Out){
		_Path = "PMNS/";
		_Filename = "generatedwith";
		_Filename += MorBStr[_B_or_M];
		_Filename += to_string(nbBitsP);
		_Filename += "delta" + to_string(_Delta);
		_Filename += "phi" + to_string(_Phi_Log);
		_Filename += "pmns.py";
	}

	/*
		if n in unbounded it will grow until a pmns is found
		if it is bounded then P will be changed until one can be found in the range
	*/
	while(true){

		
		if (_P_ORDER != -1){
			P = RandomPrime_ZZ(1 << _P_ORDER);
		}
		if (_P_Size != -1){
			P = RandomPrime_ZZ(_P_Size);
		}

		ZZ_p::init(P);

		

		//set n range
		int n_lower, n_upper;

		n_upper = _N_RANGEMAX + 1;

		if (_N != -1){
			n_lower = _N;
			n_upper = _N + 1;
		}
		
		else if (_N_RANGEMIN != -1){
			n_lower = _N_RANGEMIN;
		}
		else{
			//theoretical minimum possible n
			n_lower = (NumBits(P) / _Phi_Log)+1; //n starts at (log_2(p) / phi)) + 1
			if (_FixedRho){
				RR RHO_RR = conv<RR>(ZZ(1) << 64); 
				while(RHO_RR < rhoBound(n_lower, 2*n_lower - 1, nthRoot(P, n_lower))){
					n_lower++;
				}
			}
			else{
				while(PHI_RR < phiBound(n_lower, 2*n_lower - 1, nthRoot(P, n_lower))){
					n_lower++;
				}
			}
		}

		

		if (_Verbose){

			
			cout << "Setting P : " << P << endl;
		

			cout << "N range : [" << n_lower << " , ";
			if (n_upper > 0){
				cout << n_upper - 1;
			}
			else{
				cout << "infinity";
			}

			cout << "]" << endl;
		}
		
		for (int n = n_lower; n != n_upper; n++){

			if (_Poly_Type  == PolyType::ANY){
				_Poly_Type = PolyType::Xn_Lam;
				if (generateFromN(P, n)){
					
					return;
				}
				_Poly_Type = PolyType::Xn_X_1;
				if (generateFromN(P, n)){
					
					return;
				}
				_Poly_Type = PolyType::Xn_X2_1;
				if (generateFromN(P, n)){
					
					return;
				}
				_Poly_Type = PolyType::ANY;
			}
			else{
				if (generateFromN(P, n)){
					return;
				}
			}

			
		}
	}


}

