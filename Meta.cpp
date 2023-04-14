#include "Meta.h"

using namespace NTL;
using namespace std;

int _Phi_Log = 64;

ZZ PHI = NTL::ZZ(1) << _Phi_Log;
RR PHI_RR = NTL::conv<NTL::RR>(PHI);

//

bool _Verbose = false;
bool _Show = false;
bool _Out = false;
string _Path = "";
string _Filename = "";
int _Nb_Thread = 1;
//

NTL::ZZ _P = ZZ(-1);
int _P_ORDER = -1;
int _P_Size = -1;

int _N = -1;
int _N_RANGEMIN = -1;
int _N_RANGEMAX = -1;
int _Delta = 0;

//

PolyType _Poly_Type = PolyType::Xn_Lam;
MorB _B_or_M = MorB::M;
bool _DisableFastRootCheck = false;
bool _ExactMaxBound = false;
bool _AugmentedMinN = false;
int _FactorLimit = 0;
bool _BKZ = false;

//

bool _Xn_Lam_Lam_Set = false;

int _Xn_Lam_Lam = -1;

int _Xn_Lam_Lam_RangeMin = -1;

int _Xn_Lam_Lam_RangeMax = -1;