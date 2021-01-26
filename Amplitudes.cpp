#include <iostream> 
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <complex>
#include <ctime>

using namespace std;

const double Mp(0.93827), Mn(0.93957), Mpip(0.13957), Mpiz(0.13498);

double fRand(const double& fMin, const double& fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double F1p_f(const double& Q2)
{
	double F1p, Ge, Gm, tau;
	
	tau = Q2/(4*Mp*Mp);
	Ge = 1/((1 + Q2/0.71)*(1 + Q2/0.71));	
	Gm = 2.79/((1 + Q2/0.71)*(1 + Q2/0.71));
	F1p = (Ge + tau*Gm)/(1 + tau);	
	
	return F1p;
}

double F2p_f(const double& Q2)
{
	double F2p, Ge, Gm, tau;
	
	tau = Q2/(4*Mp*Mp);
	Ge = 1/((1 + Q2/0.71)*(1 + Q2/0.71));	
	Gm = 2.79/((1 + Q2/0.71)*(1 + Q2/0.71));
	F2p = (Gm - Ge)/(1.79*(1 + tau));	
	
	return F2p;
}

double F1n_f(const double& Q2)
{
	double F1n, Ge, Gm, tau;
	
	tau = Q2/(4*Mn*Mn);	
	Gm = -1.91/((1 + Q2/0.71)*(1 + Q2/0.71));
	Ge = -Gm*tau/(1 + 5.6*tau);
	F1n = (Ge + tau*Gm)/(1 + tau);		
	
	return F1n;
}

double F2n_f(const double& Q2)
{
	double F2n, Ge, Gm, tau;
	
	tau = Q2/(4*Mn*Mn);	
	Gm = -1.91/((1 + Q2/0.71)*(1 + Q2/0.71));
	Ge = -Gm*tau/(1 + 5.6*tau);
	F2n = (Gm - Ge)/(-1.91*(1 + tau));	
	
	return F2n;
}

double Fpi_f(const double& Q2)
{
	double Fpi;
	
	Fpi = 0.85/(1 + 0.431*Q2/(6*0.04)) + 0.15/((1 + 0.411*Q2/(12*0.04))*(1 + 0.411*Q2/(12*0.04)));	
	
	return Fpi;
}

vector<complex<double>> U_functions_n(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V2;
	double Re, Im, F1n, F1p, F2n, F2p, Fpi, s, t, u, kq, qp, nu;
	complex<double> Value, U6_1, U8_1, U5_1, U6_2, U8_2, U5_2, U6, U8;
	
	double E1, E2, k0, q0, x;	x = 13.5*sqrt(double(2)/double(137));
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mn*Mn - Mpip*Mpip)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	F1n = F1n_f(Q2); F1p = F1p_f(Q2); F2n = F2n_f(Q2); F2p = F2p_f(Q2); Fpi = Fpi_f(Q2);
	
	kq = k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta);
	s = W*W; 
	t = Mpip*Mpip - Q2 - 2*kq;
	u = Mp*Mp + Mn*Mn - Q2 + Mpip*Mpip - s - t;
	
	nu = (W*W + Q2 - Mp*Mp)/(2*Mp);
	
	qp = 2*nu*Mp - Q2 - kq;
	
	Re = -x*F1n/(u - Mp*Mp) + x*F1p/(s - Mp*Mp); //U1
	Im = -x*F2n/(u - Mp*Mp) - x*F2p/(s - Mp*Mp);	
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Value = complex<double>(0,0); //U2
	V2.push_back(Value);
	
	Re = x*Fpi*(-1/Q2 - 1/(t - Mpip*Mpip)) - x*F1n/Q2 + x*F1p/Q2; //U4
	Im = x*F2n/(u - Mp*Mp) + x*F2p/(s - Mp*Mp);
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Re = x*2*Mp*F1p/(s - Mp*Mp);
	Im = x*Q2*F2p/((s - Mp*Mp)*2*Mp);
	U5_1 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2p/(2*Mp*(s - Mp*Mp));
	U8_1 = complex<double>(Re,Im); 
	
	U6_1 = -(U5_1 - Q2*U8_1)/qp;
	
	Re = -x*2*Mp*F1n/(u - Mp*Mp);
	Im = x*Q2*F2n/((u - Mp*Mp)*2*Mp);
	U5_2 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2n/(2*Mp*(u - Mp*Mp));
	U8_2 = complex<double>(Re,Im);
	
	U6_2 = -(U5_2 - Q2*U8_2)/qp;
	
	U6 = U6_1 + U6_2;
	U8 = U8_1 + U8_2;
	
	V2.push_back(U6);
	
	Value = complex<double>(0,0); //U7
	V2.push_back(Value);
	
	V2.push_back(U8); 
	
	return V2;
}

vector<complex<double>> U_functions_p(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V2;
	double Re, Im, F1n, F1p, F2n, F2p, Fpi, s, t, u, kq, qp, nu;
	complex<double> Value, U6_1, U8_1, U5_1, U6_2, U8_2, U5_2, U6, U8;
	
	double E1, E2, k0, q0, x;	x = 13.5*sqrt(double(1)/double(137));
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mp*Mp - Mpiz*Mpiz)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	nu = (W*W + Q2 - Mp*Mp)/(2*Mp);
	
	F1n = F1n_f(Q2); F1p = F1p_f(Q2); F2n = F2n_f(Q2); F2p = F2p_f(Q2); //Fpi = Fpi_f(Q2);
	
	kq = k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta);
	s = W*W; 
	t = Mpiz*Mpiz - Q2 - 2*kq;
	u = Mp*Mp + Mn*Mn - Q2 + Mpiz*Mpiz - s - t;
	
	qp = 2*nu*Mp - Q2 - kq;
	
	Re = -x*F1p/(u - Mp*Mp) + x*F1p/(s - Mp*Mp); //U1
	Im = -x*F2p/(u - Mp*Mp) - x*F2p/(s - Mp*Mp);	
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Value = complex<double>(0,0); //U2
	V2.push_back(Value);
	
	Re = - x*F1p/Q2 + x*F1p/Q2;
	Im = x*F2p/(u - Mp*Mp) + x*F2p/(s - Mp*Mp); //U4
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Re = x*2*Mp*F1p/(s - Mp*Mp);
	Im = x*Q2*F2p/((s - Mp*Mp)*2*Mp);
	U5_1 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2p/(2*Mp*(s - Mp*Mp));
	U8_1 = complex<double>(Re,Im);
	
	U6_1 = -(U5_1 - Q2*U8_1)/qp;
	
	Re = -x*2*Mp*F1p/(u - Mp*Mp);
	Im = x*Q2*F2p/((u - Mp*Mp)*2*Mp);
	U5_2 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2p/(2*Mp*(u - Mp*Mp));
	U8_2 = complex<double>(Re,Im);
	
	U6_2 = -(U5_2 - Q2*U8_2)/qp;
	
	U6 = U6_1 + U6_2;
	U8 = U8_1 + U8_2;
	
	V2.push_back(U6);
	
	Value = complex<double>(0,0); //U7
	V2.push_back(Value);
	
	V2.push_back(U8);			
	
	return V2;
}

vector<complex<double>> A_functions_n(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V1, U;
	
	complex<double> A;
	
	double E1, E2, k0, q0, kq;
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mn*Mn - Mpip*Mpip)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	kq = k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta);
	
	U = U_functions_n(W, Q2, theta);
	
	A = U[0] - 2*Mp*U[3];
	V1.push_back(A);
	A = 0; 
	V1.push_back(A);
	A = -U[4];
	V1.push_back(A);
	A = -U[3];
	V1.push_back(A);
	A = 0;
	V1.push_back(A);
	A = U[5];
	V1.push_back(A); 
	
	return V1;
}

vector<complex<double>> A_functions_p(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V1, U;
	
	complex<double> A;
	
	double E1, E2, k0, q0, kq;
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mp*Mp - Mpiz*Mpiz)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	kq = k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta);
	
	U = U_functions_p(W, Q2, theta);
	
	A = U[0] - 2*Mp*U[3];
	V1.push_back(A);
	A = 0;
	V1.push_back(A);
	A = -U[4]; //
	V1.push_back(A);
	A = -U[3]; //
	V1.push_back(A);
	A = 0;
	V1.push_back(A);
	A = U[5];
	V1.push_back(A); 
	
	return V1;
}


vector<complex<double>> Amplitudes_Fn(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V, A;
	complex<double> F;
	
	double E1, E2, k0, q0;
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mn*Mn - Mpip*Mpip)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	A = A_functions_n(W, Q2, theta);
	
	F = sqrt((E1 + Mp)*(E2 + Mp))*(-(W - Mp)*A[0] - (k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta))*(A[2] - A[3]) - (W - Mp)*(W - Mp)*A[3] - Q2*A[5])/(8*M_PI*W);
	V.push_back(F);
	F = sqrt((E1 - Mp)*(E2 - Mp))*((W + Mp)*A[0] - (k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta))*(A[2] - A[3]) - (W + Mp)*(W + Mp)*A[3] - Q2*A[5])/(8*M_PI*W);
	V.push_back(F);
	F = (E2 + Mp)*sqrt((E1 - Mp)*(E2 - Mp))*((W*W - Mp*Mp)*A[1] - (W + Mp)*(A[2] - A[3]) - Q2*A[4])/(8*M_PI*W);
	V.push_back(F);
	F = (E2 - Mp)*sqrt((E1 + Mp)*(E2 + Mp))*(-(W*W - Mp*Mp)*A[1] - (W - Mp)*(A[2] - A[3]) + Q2*A[4])/(8*M_PI*W);
	V.push_back(F);
	F = (E1 - Mp)*sqrt((E1 + Mp)*(E2 + Mp))*(A[0] + (k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta))*(A[1] - A[4]) + (W - Mp)*A[3] + (W + Mp)*A[5])/(8*M_PI*W);
	V.push_back(F);
	F = (E1 + Mp)*sqrt((E1 - Mp)*(E2 - Mp))*(-A[0] - (k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta))*(A[1] - A[4]) + (W + Mp)*A[3] + (W - Mp)*A[5])/(8*M_PI*W);
	V.push_back(F);

	return V; 
}

vector<complex<double>> Amplitudes_Fp(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V, A;
	complex<double> F;
	
	double E1, E2, k0, q0;
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mp*Mp - Mpiz*Mpiz)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	A = A_functions_p(W, Q2, theta);
	
	F = sqrt((E1 + Mp)*(E2 + Mp))*(-(W - Mp)*A[0] - (k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta))*(A[2] - A[3]) - (W - Mp)*(W - Mp)*A[3] - Q2*A[5])/(8*M_PI*W);
	V.push_back(F);
	F = sqrt((E1 - Mp)*(E2 - Mp))*((W + Mp)*A[0] - (k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta))*(A[2] - A[3]) - (W + Mp)*(W + Mp)*A[3] - Q2*A[5])/(8*M_PI*W);
	V.push_back(F);
	F = (E2 + Mp)*sqrt((E1 - Mp)*(E2 - Mp))*((W*W - Mp*Mp)*A[1] - (W + Mp)*(A[2] - A[3]) - Q2*A[4])/(8*M_PI*W);
	V.push_back(F);
	F = (E2 - Mp)*sqrt((E1 + Mp)*(E2 + Mp))*(-(W*W - Mp*Mp)*A[1] - (W - Mp)*(A[2] - A[3]) + Q2*A[4])/(8*M_PI*W);
	V.push_back(F);
	F = (E1 - Mp)*sqrt((E1 + Mp)*(E2 + Mp))*(A[0] + (k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta))*(A[1] - A[4]) + (W - Mp)*A[3] + (W + Mp)*A[5])/(8*M_PI*W);
	V.push_back(F);
	F = (E1 + Mp)*sqrt((E1 - Mp)*(E2 - Mp))*(-A[0] - (k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta))*(A[1] - A[4]) + (W + Mp)*A[3] + (W - Mp)*A[5])/(8*M_PI*W);
	V.push_back(F);

	return V;
}

vector<double> Amplitudes_A(const double& Q2)
{
	vector<double> V;
	double Value; 
	
	Value = (21.64-12.41*Q2+1.909*Q2*Q2)/(1.0-0.4583*Q2+0.1422*Q2*Q2-0.0525*Q2*Q2*Q2+0.00931*Q2*Q2*Q2*Q2);
	V.push_back(Value/1000);
	Value = -178.45/((1.0+Q2)*(1.0+0.3457*Q2*Q2-0.087*Q2*Q2*Q2+0.00806*Q2*Q2*Q2*Q2));
	V.push_back(Value/1000);
	Value = -339.06/((1.0+Q2)*(1.0+0.3481*Q2*Q2-0.0854*Q2*Q2*Q2+0.00758*Q2*Q2*Q2*Q2));
	V.push_back(Value/1000);
	
	return V;
}

double Sections_from_CGLN_p(vector<complex<double>>& V, const double& W, const double& Q2, const int& phi, const double& theta, const double& E0)
{
	double Section, Rt, Rl, Rtl, Rtt, Rtl2, L(19.732688), nu, CC, C, eps, Gamma_flux;
	complex<double> F1, F2, F3, F4, F5, F6;
	
	double Epi = (W*W + Mpiz*Mpiz - Mp*Mp)/(2*W);
	double Ppi = sqrt(Epi*Epi - Mpiz*Mpiz);
	
	int polarization(0);
	
	nu = (W*W + Q2 - Mp*Mp)/(2*Mp);
	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));	
	CC = sqrt(Q2)/nu;
	
	C = 2*W*Ppi/sqrt((pow((W*W - Mp*Mp), 2) + Q2*(2*Mp*Mp + Q2)));
	
	F1 = V[0];
	F2 = V[1];
	F3 = V[2];
	F4 = V[3];
	F5 = V[4];
	F6 = V[5];

	Rt = pow(L, 2)*(abs(F1)*abs(F1) + abs(F2)*abs(F2) + 0.5*pow(sin(theta), 2)*(abs(F3)*abs(F3) + abs(F4)*abs(F4)) - real(2*cos(theta)*conj(F1)*F2 - pow(sin(theta),2)*(conj(F1)*F4 + conj(F2)*F3 + cos(theta)*conj(F3)*F4)));
	Rl = CC*CC*pow(L, 2)*(abs(F5)*abs(F5) + abs(F6)*abs(F6) + 2*cos(theta)*real(conj(F5)*F6));
	Rtl = CC*pow(L, 2)*(-sin(theta)*real((conj(F2) + conj(F3) + cos(theta)*conj(F4))*F5 + (conj(F1) + conj(F4) + cos(theta)*conj(F3))*F6));
	Rtt = pow(L, 2)*pow(sin(theta), 2)*(0.5*(abs(F3)*abs(F3) + abs(F4)*abs(F4)) + real(conj(F1)*F4 + conj(F2)*F3 + cos(theta)*conj(F3)*F4));
	Rtl2 = CC*pow(L, 2)*(-sin(theta)*imag((conj(F2) + conj(F3) + cos(theta)*conj(F4))*F5 + (conj(F1) + conj(F4) + cos(theta)*conj(F3))*F6));

	Section = C*Rt + eps*C*Rl + sqrt(2*eps*(1 + eps))*Rtl*cos(phi)*C + eps*C*Rtt*cos(2*phi) + polarization*C*sqrt(2*eps*(1 - eps))*Rtl2*sin(phi);
	
	cout << "Born/Proton/Check (mubn/sr)" << endl;
	
	cout << "\nS_t = " << 4*M_PI*C*Rt << endl;
	cout << "S_l = " << 4*M_PI*C*Rl << endl;
	cout << "S_tl = " << 4*M_PI*C*Rtl << endl;
	cout << "S_tt = " << 4*M_PI*C*Rtt << endl;
	cout << "S_tl2 = " << 4*M_PI*C*Rtl2 << "\n\n" << endl;
	
	
	Gamma_flux = W*(W*W - Mp*Mp)/(137*4*M_PI*Mp*Mp*E0*E0*(1 - eps)*Q2);

	Section = Gamma_flux*Section;

	return Section;
}

double Sections_from_CGLN_n(vector<complex<double>>& V, const double& W, const double& Q2, const int& phi, const double& theta, const double& E0)
{
	double Section, Rt, Rl, Rtl, Rtt, Rtl2, L(19.732688), nu, CC, C, eps, Gamma_flux;
	complex<double> F1, F2, F3, F4, F5, F6;
	
	double Epi = (W*W + Mpip*Mpip - Mp*Mp)/(2*W);
	double Ppi = sqrt(Epi*Epi - Mpip*Mpip);
	
	int polarization(0);
	
	nu = (W*W + Q2 - Mp*Mp)/(2*Mp);
	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));	
	CC = sqrt(Q2)/nu;
	
	C = 2*W*Ppi/sqrt((pow((W*W - Mp*Mp), 2) + Q2*(2*Mp*Mp + Q2)));
	
	F1 = V[0];
	F2 = V[1];
	F3 = V[2];
	F4 = V[3];
	F5 = V[4];
	F6 = V[5];

	Rt = pow(L, 2)*(abs(F1)*abs(F1) + abs(F2)*abs(F2) + 0.5*pow(sin(theta), 2)*(abs(F3)*abs(F3) + abs(F4)*abs(F4)) - real(2*cos(theta)*conj(F1)*F2 - pow(sin(theta),2)*(conj(F1)*F4 + conj(F2)*F3 + cos(theta)*conj(F3)*F4)));
	Rl = CC*CC*pow(L, 2)*(abs(F5)*abs(F5) + abs(F6)*abs(F6) + 2*cos(theta)*real(conj(F5)*F6));
	Rtl = CC*pow(L, 2)*(-sin(theta)*real((conj(F2) + conj(F3) + cos(theta)*conj(F4))*F5 + (conj(F1) + conj(F4) + cos(theta)*conj(F3))*F6));
	Rtt = pow(L, 2)*pow(sin(theta), 2)*(0.5*(abs(F3)*abs(F3) + abs(F4)*abs(F4)) + real(conj(F1)*F4 + conj(F2)*F3 + cos(theta)*conj(F3)*F4));
	Rtl2 = CC*pow(L, 2)*(-sin(theta)*imag((conj(F2) + conj(F3) + cos(theta)*conj(F4))*F5 + (conj(F1) + conj(F4) + cos(theta)*conj(F3))*F6));

	cout << "\n\nBorn/Neutron/Check (mubn/sr)\n" << endl;
	
	cout << "S_t = " << 4*M_PI*C*Rt << endl;
	cout << "S_l = " << 4*M_PI*C*Rl << endl;
	cout << "S_tl = " << 4*M_PI*C*Rtl << endl;
	cout << "S_tt = " << 4*M_PI*C*Rtt << endl;
	cout << "S_tl2 = " << 4*M_PI*C*Rtl2 << "\n\n" << endl;

	Section = C*Rt + eps*C*Rl + sqrt(2*eps*(1 + eps))*Rtl*cos(phi)*C + eps*C*Rtt*cos(2*phi) + polarization*C*sqrt(2*eps*(1 - eps))*Rtl2*sin(phi);
	
	Gamma_flux = W*(W*W - Mp*Mp)/(137*4*M_PI*Mp*Mp*E0*E0*(1 - eps)*Q2);

	Section = Gamma_flux*Section;

	return Section;
}

int main(int argc, char **argv) // W, Q2, cos(theta) generation
{
	srand(time(NULL));
	
	vector<complex<double>> Amp1, Amp3;
	vector<double> Amp2;
	
	double S_t, S_l, S, eps, Gamma_flux, S1, S2;	
	double W, Q2, theta, E0(6.5), nu, phi; 
	
	W = fRand(1.2, 1.3); 
	Q2 = fRand(0.3, 1); 
	theta = fRand(0, M_PI); 
	phi = fRand(0, 2*M_PI);
	
	Amp1 = Amplitudes_Fn(W, Q2, theta);
	Amp3 = Amplitudes_Fp(W, Q2, theta);
	Amp2 = Amplitudes_A(Q2);

	nu =  (W*W + Q2 - Mp*Mp)/(2*Mp);
	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));
	Gamma_flux = W*(W*W - Mp*Mp)/(137*4*M_PI*Mp*Mp*E0*E0*(1 - eps)*Q2);
	
	S_t = 19.732688*19.732688*2*Mp*(abs(Amp2[1])*abs(Amp2[1]) + abs(Amp2[2])*abs(Amp2[2]))/(1.232*0.116);
	S_l = 19.732688*19.732688*4*Mp*(abs(Amp2[0])*abs(Amp2[0]))/(1.232*0.116);
	S = Gamma_flux*(S_t + eps*S_l);
	
	cout << "\n\nDelta(1232)/Check (mubn/sr)\n" << endl;
	
	cout << "S_t_p = " << 8*M_PI*S_t/5 << endl;
	cout << "S_l_p = " << 8*M_PI*S_l/5 << endl;
	cout << "S_t_n = " << 12*M_PI*S_t/5 << endl;
	cout << "S_l_n = " << 12*M_PI*S_l/5 << endl;
	
	S1 = Sections_from_CGLN_n(Amp1, W, Q2, phi, theta, E0);
	S2 = Sections_from_CGLN_p(Amp3, W, Q2, phi, theta, E0);
	
	cout << "-----------------------------------------------" << endl;
	cout << "\tKinematic values\n\nW\t\t=\t" << W << "\t\tGeV\nQ2\t\t=\t" << Q2 << "\t\tGeV^2\nCos(Theta)\t=\t" << cos(theta) << endl;
	cout << "phi\t\t=\t" << phi*180/M_PI << "\nE_0\t\t=\t" << E0 << "\t\tGeV" << endl; 
	cout << "-----------------------------------------------" << endl;
	cout << "\tAmplitudes\n" << endl;
	cout << "(1) CGLN Amplitudes for Pi^{+}n final system (GeV^-1):\n" << endl;
	cout << "\tF1\t=\t" << Amp1[0].real() << "\t+\t" << Amp1[0].imag() << endl;	
	cout << "\tF2\t=\t" << Amp1[1].real() << "\t+\t" << Amp1[1].imag() << endl;
	cout << "\tF3\t=\t" << Amp1[2].real() << "\t+\t" << Amp1[2].imag() << endl;
	cout << "\tF4\t=\t" << Amp1[3].real() << "\t+\t" << Amp1[3].imag() << endl;
	cout << "\tF5\t=\t" << Amp1[4].real() << "\t+\t" << Amp1[4].imag() << endl;
	cout << "\tF6\t=\t" << Amp1[5].real() << "\t+\t" << Amp1[5].imag() << endl;
	cout << "\n(2) CGLN Amplitudes for Pi^{0}p final system (GeV^-1):\n" << endl;
	cout << "\tF1\t=\t" << Amp3[0].real() << "\t+\t" << Amp3[0].imag() << endl;	
	cout << "\tF2\t=\t" << Amp3[1].real() << "\t+\t" << Amp3[1].imag() << endl;
	cout << "\tF3\t=\t" << Amp3[2].real() << "\t+\t" << Amp3[2].imag() << endl;
	cout << "\tF4\t=\t" << Amp3[3].real() << "\t+\t" << Amp3[3].imag() << endl;
	cout << "\tF5\t=\t" << Amp3[4].real() << "\t+\t" << Amp3[4].imag() << endl;
	cout << "\tF6\t=\t" << Amp3[5].real() << "\t+\t" << Amp3[5].imag() << endl;
	cout << "\n(3) Helicity Amplitudes for Delta(1232) resonance:\n" << endl;
	cout << "\tS12\t=\t" << Amp2[0] << endl;	
	cout << "\tA12\t=\t" << Amp2[1] << endl;
	cout << "\tA32\t=\t" << Amp2[2] << endl;
	cout << "\n(4) Diff. cross-setions:\n" << endl; 
	cout << "\tdS_(Delta_1232)\t=\t" << S << "\t\tmubn" << endl;	
	cout << "\tdS_(Pi^{0}p)\t=\t" << S2 << "\tmubn" << endl;
	cout << "\tdS_(Pi^{+}n)\t=\t" << S1 << "\tmubn" << endl;

	
	Amp1.clear(); Amp2.clear(); Amp3.clear();

	return 0;
}
