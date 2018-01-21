#include "INT_quadratureIntegration.h"
//*********************************************************
//*********************************************************
//*********************************************************
//			INTEGRATE
// Integrates input function f(i) wrt r, from l to m. This program contains wronskian drdt
//Uses an nquad-point quadrature formula.. for nquad=1->14 (any integer)
//double integrate(double *f, int l, int m)
double INT_integrate(std::vector<double> f, std::vector<double> w, int l,
  int m, int nquad)
/*

f(x) dx -> f(u(t)) w(t) dt

w = du/dt = dr/dt is the wronskian

XXX overload so w is optional!
XXX overload so can use floats?

*/
{

//Note: normally, NGP is so large, that there is no difference between nquad=1 to 14!


//re-arange limits to make less fiddly
// so I don't need to ..
	int negInt=1;
	if (l>m){
		int temp=l;
		l=m;
		m=temp;
		negInt=-1;
	}
	if (m>=NGP){m=NGP-1;}
	if (l<0){l=0;}


// Defines the `nquad'-point quadrature integration coeficents.
// nquad can take any integer from 1 to 14 (=1 implies trapazoid rule)
	double ic[nquad];
	if (nquad==1){
		double c[1]={1};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==2){
		double c[2]={1,2};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==3){
		double c[3]={9,28,23};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==4){
		double c[4]={8,31,20,25};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==5){
		double c[5]={475,1902,1104,1586,1413};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==6){
		double c[6]={459,1982,944,1746,1333,1456};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==7){
		double c[7]={36799,176648,54851,177984,89437,130936,119585};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==8){
		double c[8]={35584,185153,29336,220509,46912,156451,111080,122175};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==9){
		double c[9]={2082753,11532470,261166,16263486,-1020160,12489922,5095890,
		7783754,7200319};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==10){
		double c[10]={2034625,11965622,-1471442,20306238,-7084288,18554050,1053138,
		9516362,6767167,7305728};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==11){
		double c[11]={262747265,1637546484,-454944189,3373884696,-2145575886,3897945600,
		-1065220914,1942518504,636547389,1021256716,952327935};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==12){
		double c[12]={257696640,1693103359,-732728564,4207237821,-3812282136,6231334350,
		-3398609664,3609224754,-196805736,1299041091,896771060,963053825};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==13){
		double c[13]={1382741929621,9535909891802,-5605325192308,28323664941310,
		-32865015189975,53315213499588,-41078125154304,39022895874876,-13155015007785,
		12465244770050,3283609164916,5551687979302,5206230892907};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else if (nquad==14){
		double c[14]={1360737653653,9821965479386,-7321658717812,34616887868158,
		-48598072507095,81634716670404,-78837462715392,76782233435964,-41474518178601,
		28198302087170,-3009613761932,7268021504806,4920175305323,5252701747968};
		for (int i=0;i<nquad;i++){ic[i]=c[i];}
	}
	else {
		printf("FAILURE: Wrong order for integration.. check nquad\n");
		return 0;
	}

	double dd[14]={2,2,24,24,1440,1440,120960,120960,7257600,
					7257600,958003200,958003200,5230697472000,5230697472000};


	// wronskian
  //XXX make optional!
	for (int i=l; i<=m; i++){
		f[i]=f[i]*w[i]; //XXX double check - is this OK??
	}

	// Defines/calculates the integral
	double a=0;
	int ll=l;
	int mm=m;
	for (int i=0; i<nquad; i++){
		a=a+ic[i]*(f[ll]+f[mm]);
		ll=ll+1;
		mm=mm-1;
	}
	a=a/dd[nquad-1];
	for (int i=l; i<=m; i++){
		a=a+f[i];
	}
	a=a*h*negInt;

	// // rectangle method (for tests)
	// double b=0;
	// for (int i=l; i<m; i++){
	// 	b=b+f[i];
	// }
	// b=b*h;

	return a;
}		// END integrate
