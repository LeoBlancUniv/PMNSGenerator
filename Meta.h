#ifndef META_H
#define META_H

#include <NTL/ZZ.h>

#include <string.h>

#include "Poly.h"
#include "FindM.h"

/*
	file to hold global variable and define
*/

#define FACTOR_LIMIT_DEFAULT 80
#define LAM_START_DEFAULT 1



enum PolyType : int;
enum MorB : int;

//

extern int _Phi_Log;

extern NTL::ZZ PHI;
extern NTL::RR PHI_RR;

//

extern bool _Verbose;
extern bool _SVerbose;
extern bool _Show;
extern bool _Out;
extern std::string _Path;
extern std::string _Filename;
extern int _Nb_Thread;
//

extern NTL::ZZ _P;
extern int _P_ORDER;
extern int _P_Size;


extern int _N;
extern int _N_RANGEMIN;
extern int _N_RANGEMAX;
extern int _Delta;

//

extern PolyType _Poly_Type;
extern MorB _B_or_M;
extern bool _DisableFastRootCheck;
extern bool _ExactMaxBound;
extern bool _AugmentedMinN;
extern bool _ForceIrreducible;
extern int _FactorLimit;
extern bool _BKZ;
extern bool _FixedRho;


//

extern bool _Xn_Lam_Lam_Set;

extern int _Xn_Lam_Lam;

extern int _Xn_Lam_Lam_RangeMin;

extern int _Xn_Lam_Lam_RangeMax;


#endif