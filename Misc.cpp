#include "Misc.h"

using namespace std;
using namespace NTL;
using namespace fplll;








Vec<int> primeFactor(int n){
	/*
		factor <n> 
		super basic way,  do not use on big number
	*/
	Vec<int> Factor_vec;

	for (int i = 2; i < n+1; i++){
		while(n%i == 0){
			Factor_vec.append(i);
			n /= i;
		}
	}

	return Factor_vec;

}

void convertZZtoMPZ(mpz_t Dest, const ZZ& Src){
	/*
		convert a NTL ZZ number to a GMP mpz one 
	*/
	
	const ZZ_limb_t* data = ZZ_limbs_get(Src);
	
	int size = Src.size();

	mpz_init(Dest);

	mpz_import(Dest, size, -1, sizeof(data[0]), 0, 0, data);

	if (Src < 0){
		Dest->_mp_size = -Dest->_mp_size;
	}
}

void convertMPZtoZZ(ZZ& Dest, const mpz_t Src){
	/*
		convert a GMP mpz number to a NTL ZZ one 
	*/
	const mp_limb_t* data = mpz_limbs_read(Src);

	int size = mpz_size(Src);

	ZZ_limbs_set(Dest, data, size);

	if (Src->_mp_size < 0){
		Dest = -Dest;
	}
}

void convertMatZZtoMatMPZ(ZZ_mat<mpz_t>& Dest, const Mat<ZZ>& Src){
	/*
		convert a ZZ matrix to a GMP one
	*/
	int nBRowLine = Src.NumRows();

	ZZ_mat<mpz_t> Ret(nBRowLine, nBRowLine);

	for (int i = 0; i < nBRowLine; i++){
		for (int j = 0; j < nBRowLine; j++){

			mpz_t Tmp;

			convertZZtoMPZ(Tmp, Src[i][j]);

			Ret[i][j] = Tmp;
		}
	}

	Dest = Ret;
}

void convertMatMPZtoMatZZ(Mat<ZZ>& Dest, const ZZ_mat<mpz_t>& Src){
	/*
		convert a GMP matrix to a ZZ one
	*/
	int nBRowLine = Src.get_rows();

	Mat<ZZ> Ret;

	Ret.SetDims(nBRowLine, nBRowLine);

	for (int i = 0; i < nBRowLine; i++){
		for (int j = 0; j < nBRowLine; j++){

			ZZ Tmp;

			convertMPZtoZZ(Tmp, Src[i][j].get_data());

			Ret[i][j] = Tmp;
		}
	}
	Dest = Ret;
}

void freeMatMPZ(ZZ_mat<mpz_t>& Mat){
	/*
		free a GMP mpz matrix
	*/
	int nBRowLine = Mat.get_rows();

	

	for (int i = 0; i < nBRowLine; i++){
		for (int j = 0; j < nBRowLine; j++){
			
			mpz_clear(Mat[i][j].get_data());
			Mat[i][j].get_data()->_mp_d = 0;
			
		}
	}

}



void writeVec(ofstream& File, Vec<ZZ> Vec){
	/*
		write the vector <Vec> into the <File> in the formar "[a, b, c, .., z]"
	*/
	int len = Vec.length();


	File << "[";
	for (int i = 0; i < len; i++){
		File << Vec[i];

		if (i != len - 1){
			File << ", ";
		}
	}

	File << "]";
}

void writeMatrix(ofstream& File, Mat<ZZ> Mat){
	/*
		write the matrix <Mat> into the <File> in the format "[vec0, vec1, ..., vecn]"
	*/
	int len = Mat.NumRows();

		File << "[";

		for (int i = 0; i < len; i++){

			Vec<ZZ> TmpVec = Mat[i];

			writeVec(File, TmpVec);

			if (i != len - 1){
				File << ", ";
			}
		}
		File << "]";
}

void writeToFile(string filename, const ZZ& P, const int n, const ZZ_p& Gamma,
				 const ZZ& Rho, const ZZX& E, const Mat<ZZ>& Base,
				 const Mat<ZZ>& Base_inv, const ZZX& M, const ZZX& M_inv){
	/*
		write a pmns into a file, format is
		add "# PHI <PHI> , Delta <Delta>" and "pmnsdict = {}" if it's a new file
		then on one line
		"pmnsdict[<P>] = (<P>, <n>, <Gamma>, relevent coef of <E>, Rho, <M> or <B>, <M_inv> or <B_inv>)"
	*/
	bool exist = true;

	if (FILE *file = fopen(filename.c_str(), "r")) { // causes an unused variable error
		fclose(file);
	}
	else{
		exist = false;
	}

	ofstream File;

	File.open(filename, ios::app);

	if (!exist){
		File << "# PHI : " << _Phi_Log << " , Delta : " << _Delta << endl;
		File << "pmnsdict = {}" << endl;
	}

	File << "pmnsdict[" << P << "] = (" 
	<< P << ", " 
	<< n << ", "
	<< Gamma << ", " ;

	Vec<ZZ> Coeff_vec;

	switch (_Poly_Type){
		case Xn_Lam : 
			File << -ConstTerm(E);
		break;

		case Xn_X_1 : 
			Coeff_vec.append(ConstTerm(E));
			Coeff_vec.append(coeff(E, 1));
		break;

		case Xn_X2_1 : 
			Coeff_vec.append(ConstTerm(E));

			for (int i = 1; i < n/2; i++){
				Coeff_vec.append(ZZ(0));
			}

			Coeff_vec.append(coeff(E, n/2));

		break;

		default : cout << "ERROR in _Poly_Type" << endl; break;
	}

	if (Coeff_vec.length() > 0){
		writeVec(File, Coeff_vec);
	}

	File << ", ";
	
	File << Rho << ", ";

	if (_B_or_M == MorB::B){
		writeMatrix(File, Base);
		File << ", ";
		writeMatrix(File, Base_inv);
		
	}
	else{
		writeVec(File, conv<Vec<ZZ>>(M));
		File << ", ";
		writeVec(File, conv<Vec<ZZ>>(M_inv));
	}

	File << ")" << endl;


	
	File.close();

}