#include <iostream>
#include <string.h>

#include <NTL/ZZ.h>


#include "MNSGenerator.h"
#include "Meta.h"
#include "Poly.h"

using namespace std;
using namespace NTL;

int main(int argc, char const *argv[])
{
	if (argc == 1){
		cerr << "too little arg, see README" << endl;
		return -1;
	}

	int i = 1;

	while(i < argc){
		string current = argv[i];
		
		//preset

		if (current == "-EnsureSmallest"){
			_Poly_Type = PolyType::ANY;
			_DisableFastRootCheck = true;
			_FactorLimit = -1;

			i++;
		}

		else if (current == "-FastHeuristic"){
			_Poly_Type = PolyType::Xn_Lam;
			_FactorLimit = 60;
			_ExactMaxBound = true;
			_AugmentedMinN = true;
			i++;
		}

		//misc option

		else if (current == "-V"){
			_Verbose = true;
			_Show = true;
			i++;
		}

		else if (current == "-SV"){
			_SVerbose = true;
			_Verbose = true;
			_Show = true;
			i++;
		}

		else if (current == "-Show"){
			_Show = true;
			i++;
		}

		else if (current == "-Out"){
			_Out = true;
			i++;
		}

		else if (current == "-Nb_Thread"){
			_Nb_Thread = stoi(argv[i+1]);
			i += 2;
		}
		
		//param option
		else if (current == "-P"){
			_P = conv<ZZ>(argv[i+1]);
			i += 2;
		}
		else if (current == "-P_order"){
			_P_ORDER = stoi(argv[i+1]);
			i += 2;
		}
		else if (current == "-P_size"){
			_P_Size = stoi(argv[i+1]);
			i += 2;
		}

		else if (current == "-N"){
			_N = stoi(argv[i+1]);
			i += 2;
		}
		else if (current == "-N_RangeMin"){
			_N_RANGEMIN = stoi(argv[i+1]);
			i += 2;
		}
		else if (current == "-N_RangeMax"){
			_N_RANGEMAX = stoi(argv[i+1]);
			i += 2;
		}

		else if (current == "-Phi"){
			_Phi_Log = stoi(argv[i+1]);
			PHI = NTL::ZZ(1) << _Phi_Log;
			PHI_RR = NTL::conv<NTL::RR>(PHI);
			i += 2;
		}

		else if (current == "-Delta"){
			_Delta = stoi(argv[i+1]);
			i += 2;
		}

		//logic option
		
		else if (current == "-Poly_Type"){
			string next = argv[i+1];
			if (next == "Xn_Lam"){
				_Poly_Type = PolyType::Xn_Lam;
			}
			else if (next == "Xn_X_1"){
				_Poly_Type = PolyType::Xn_X_1;
			}
			else if (next == "Xn_X2_1"){
				_Poly_Type = PolyType::Xn_X2_1;
			}
			else if (next == "ANY"){
				_Poly_Type = PolyType::ANY;
			}
			else{
				cerr << "param \"" << argv[i+1] << "\"" << " not recognised" << endl;
				return -1;
			}
			i += 2;
		}

		else if (current == "-B_or_M"){
			string next = argv[i+1];
			if (next == "B"){
				_B_or_M = MorB::B;
			}
			else if (next == "M"){
				_B_or_M = MorB::M;
			}
			else{
				cerr << "param \"" << argv[i+1] << "\"" << " not recognised" << endl;
				return -1;
			}
			i += 2;
		}

		else if (current == "-DisableFastRootCheck"){
			_DisableFastRootCheck = true;
			i++;
		}

		else if (current == "-ExactMaxBound"){
			_ExactMaxBound = true;
			i++;
		}

		else if (current == "-AugmentedMinN"){
			_AugmentedMinN = true;
			i++;
		}

		else if (current == "-FactorLimit"){
			int next = stoi(argv[i+1]);
			_FactorLimit = next;
			i += 2;
		}

		else if (current == "-BKZ"){
			_BKZ = true;
			i++;
		}

		else if (current == "-Xn_Lam_Lam_Set"){
			_Xn_Lam_Lam_Set = true;
			_Xn_Lam_Lam = stoi(argv[i+1]);

			i += 2;
		}

		else if (current == "-Xn_Lam_Lam_RangeMin"){
			_Xn_Lam_Lam_RangeMin = stoi(argv[i+1]);
			i += 2;
		}

		else if (current == "-Xn_Lam_Lam_RangeMax"){
			_Xn_Lam_Lam_RangeMax = stoi(argv[i+1]);
			i += 2;
		}

		else{
			cerr << "param \"" << argv[i] << "\"" << " not recognised" << endl;
			return -1;
		}
	}

	//some logic check

	



	//we have set all the meta parameter 
	MNSGenerator();

	

	return 0;
}