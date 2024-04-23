
/*
This test factorizes a given polynomial and finds its roots
*/
//*********************************************************************
#include"Rand.h"

// this function: "findroots" finds roots of a polynomial. It is a wrapper function that calls
// NTL's function: FindRoots.
// findroots returns an arrary of roots and also modifies the "number_of_roots" found.

bigint* findroots(bigint *coeff, int coeff_size, int& number_of_roots, bigint pubmoduli, ZZ_pX P){

	int counter_roots = 0;
	bigint *res;
	res = (mpz_t*)malloc(coeff_size * sizeof(mpz_t));
	ZZ one(1);
	for(int j = 0; j < coeff_size; j++){
		char * tmp = mpz_get_str(NULL, 10, coeff[j]);
		ZZ_p dd = to_ZZ_p(conv<ZZ> (tmp));
		SetCoeff(P, j, dd);
	}
	ZZ_p a = LeadCoeff(P);
	ZZ aa = rep(a);
	if(aa > one){
		MakeMonic(P);
	}
	Vec< Pair < ZZ_pX, long > > factors;
	CanZass(factors, P);
	vec_ZZ_p root;
	for(int j = 0; j < factors.length(); j++){
		if(factors[j].a.rep.length() == 2){
			root = FindRoots(factors[j].a);
			for(int k = 0; k < root.length(); k++){
				stringstream ss;
				ss << root[k];
				string tmpm = ss.str();
				char cv[tmpm.length()];
				strcpy(cv, tmpm.c_str());
				mpz_init_set_str(res[counter_roots], cv, 10);
				counter_roots++;
				tmpm.clear();
			}
		}
	}
	number_of_roots = counter_roots;
	return res;
}


int main(){
	bigint *pu_moduli;
	double start_enc;
  double end_enc;
  double diff_enc = 0;
	float count = 0;
	int number_of_roots = 0;
	Random rd;
	int loop_size = 100;// number of experiments
	int bit_size = 256;//or 128
	mpz_init_set_str(pu_moduli[0],"247742461804815734743746154725192129101",10);// this public moduli is of size 128 bits
	//mpz_init_set_str(pu_moduli[0],"189415934696602420616687117501983124694655143342576097933654834670629027463",10);// uncommnet this for 256-bit field
	int coeff_size = 11;

	// you can use the following two lines to generate a random prime number of length bit_size
	// then you can set pu_moduli[0] to this prime number for the later executions
	// mpz_nextprime(pu_moduli[0], pu_moduli[0]);
	// pu_moduli = rd.gen_randSet(1, bit_size);
	bigint *coeff;
	coeff = (mpz_t*)malloc(coeff_size * sizeof(mpz_t));
	// setting the random coefficients of polynomials. For a polynomial of degree x, it is
	// recommended to set the coeff[x] to 1, as it is closer to our schemes setting.
	// below we generate a polynomial of degree 10, as coeff[10].x_10+ coeff[9].x_9+...+coeff[0] mod (pu_moduli[0])
	// to reduce the degree of the polynomial, please comment out (or add) related coefficients.
	// also make sure variable: coeff_size is adjusted accordingly
	mpz_init_set_str(coeff[0], "6955865768787685865767656757", 10);
	mpz_init_set_str(coeff[1], "934987175705414987679878768765", 10);
	mpz_init_set_str(coeff[2], "20091969869909009898696", 10);
	mpz_init_set_str(coeff[3], "8769670916986986932869862543986", 10);
	mpz_init_set_str(coeff[4], "06429799869876435656681", 10);
	mpz_init_set_str(coeff[5], "68768897987768768768761", 10);
	mpz_init_set_str(coeff[6], "9878986796896979879698698691", 10);
	mpz_init_set_str(coeff[7], "43681519987613243234", 10);
	mpz_init_set_str(coeff[8], "097097091", 10);
	mpz_init_set_str(coeff[9], "24632908099863456643773809801", 10);
	mpz_init_set_str(coeff[10], "1", 10);

	char * tmp_mod = mpz_get_str(NULL, 10, pu_moduli[0]);
	ZZ p = to_ZZ(tmp_mod);
	bigint*rees;
	ZZ_p::init(p);
	ZZ_pX P;
	for(int i = 0; i < loop_size; i++){
		start_enc = 0;
		start_enc =clock();
 		rees = findroots( coeff, coeff_size,  number_of_roots,  pu_moduli[0], P);
		end_enc = 0;
		end_enc = clock();
		diff_enc = 0;
		diff_enc = end_enc - start_enc;
		count += diff_enc/(double) CLOCKS_PER_SEC;
}
cout<<"\n average time for "<<loop_size<<" experiments: "<<count/loop_size<<endl;
cout<<"\n pu_moduli[0]: "<<pu_moduli[0]<<endl;
cout<<"\n number_of_roots: "<<number_of_roots<<endl;
cout<<"\n rees[0]: "<<rees[0]<<endl;
return 0;
}



//**********************************************************************
