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
	double Re, Im, F1n, F1p, F2n, F2p, Fpi, s, t, u, kq;
	complex<double> Value, U7_1, U8_1, U5_1, U7_2, U8_2, U5_2, U7, U8;
	
	double E1, E2, k0, q0, x;	x = sqrt(double(2)/double(137));
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mn*Mn - Mpip*Mpip)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	F1n = F1n_f(Q2); F1p = F1p_f(Q2); F2n = F2n_f(Q2); F2p = F2p_f(Q2); Fpi = Fpi_f(Q2);
	
	kq = k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta);
	s = W*W; 
	t = Mpip*Mpip - Q2 - 2*kq;
	u = Mp*Mp + Mn*Mn - Q2 + Mpip*Mpip - s - t;
	
	Re = -x*F1n/(u - Mp*Mp) + x*F1p/(s - Mp*Mp); 
	Im = -x*F2n/(u - Mp*Mp) - x*F2p/(s - Mp*Mp);	
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Value = complex<double>(0,0);
	V2.push_back(Value);
	
	Re = x*Fpi*(-1/Q2 - 1/(t - Mpip*Mpip)) - x*F1n/Q2 + x*F1p/Q2;
	Im = x*F2n/(u - Mp*Mp) + x*F2p/(s - Mp*Mp);
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Value = complex<double>(0,0);
	V2.push_back(Value);
	
	Re = x*2*Mp*F1p/(s - Mp*Mp);
	Im = x*Q2*F2p/((s - Mp*Mp)*2*Mp);
	U5_1 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2p/(2*Mp*(s - Mp*Mp));
	U8_1 = complex<double>(Re,Im);
	
	U7_1 = -(U5_1 - Q2*U8_1)/kq;
	
	Re = -x*2*Mp*F1n/(u - Mp*Mp);
	Im = x*Q2*F2n/((u - Mp*Mp)*2*Mp);
	U5_2 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2n/(2*Mp*(u - Mp*Mp));
	U8_2 = complex<double>(Re,Im);
	
	U7_2 = -(U5_2 - Q2*U8_2)/kq;
	
	U7 = U7_1 + U7_2;
	U8 = U8_1 + U8_2;
	
	V2.push_back(U7);
	V2.push_back(U8); 
	
	return V2;
}

vector<complex<double>> U_functions_p(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V2;
	double Re, Im, F1n, F1p, F2n, F2p, Fpi, s, t, u, kq;
	complex<double> Value, U7_1, U8_1, U5_1, U7_2, U8_2, U5_2, U7, U8;
	
	double E1, E2, k0, q0, x;	x = sqrt(double(1)/double(137));
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mp*Mp - Mpiz*Mpiz)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	F1n = F1n_f(Q2); F1p = F1p_f(Q2); F2n = F2n_f(Q2); F2p = F2p_f(Q2); Fpi = Fpi_f(Q2);
	
	kq = k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta);
	s = W*W; 
	t = Mpiz*Mpiz - Q2 - 2*kq;
	u = Mp*Mp + Mn*Mn - Q2 + Mpiz*Mpiz - s - t;
	
	Re = -x*F1p/(u - Mp*Mp) + x*F1p/(s - Mp*Mp);
	Im = -x*F2p/(u - Mp*Mp) - x*F2p/(s - Mp*Mp);	
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Value = complex<double>(0,0);
	V2.push_back(Value);
	
	Re = - x*F1p/Q2 + x*F1p/Q2;
	Im = x*F2p/(u - Mp*Mp) + x*F2p/(s - Mp*Mp);
	Value = complex<double>(Re,Im);
	V2.push_back(Value);
	
	Value = complex<double>(0,0);
	V2.push_back(Value);
	
	Re = x*2*Mp*F1p/(s - Mp*Mp);
	Im = x*Q2*F2p/((s - Mp*Mp)*2*Mp);
	U5_1 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2p/(2*Mp*(s - Mp*Mp));
	U8_1 = complex<double>(Re,Im);
	
	U7_1 = -(U5_1 - Q2*U8_1)/kq;
	
	Re = -x*2*Mp*F1p/(u - Mp*Mp);
	Im = x*Q2*F2p/((u - Mp*Mp)*2*Mp);
	U5_2 = complex<double>(Re,Im);
	
	Re = 0;
	Im = x*F2p/(2*Mp*(u - Mp*Mp));
	U8_2 = complex<double>(Re,Im);
	
	U7_2 = -(U5_2 - Q2*U8_2)/kq;
	
	U7 = U7_1 + U7_2;
	U8 = U8_1 + U8_2;
	
	V2.push_back(U7);
	V2.push_back(U8);			
	
	return V2;
}

vector<complex<double>> A_functions_n(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V1, U;
	
	complex<double> A;
	
	double E1, E2, k0, q0;
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mn*Mn - Mpip*Mpip)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	U = U_functions_n(W, Q2, theta);
	
	A = U[0] - 2*Mp*U[3];
	V1.push_back(A);
	A = U[1]/(k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta));
	V1.push_back(A);
	A = -U[4];
	V1.push_back(A);
	A = -U[3];
	V1.push_back(A);
	A = (U[0] + U[2])/(k0*q0 - sqrt(k0*k0 - Mpip*Mpip)*sqrt(q0*q0 + Q2)*cos(theta));
	V1.push_back(A);
	A = U[5];
	V1.push_back(A); 
	
	return V1;
}

vector<complex<double>> A_functions_p(const double& W, const double& Q2, const double& theta)
{
	vector<complex<double>> V1, U;
	
	complex<double> A;
	
	double E1, E2, k0, q0;
	
	E1 = (W*W + Mp*Mp + Q2)/(2*W);
	E2 = (W*W + Mp*Mp - Mpiz*Mpiz)/(2*W);
	k0 = W - E2;
	q0 = W - E1;	
	
	U = U_functions_p(W, Q2, theta);
	
	A = U[0] - 2*Mp*U[3];
	V1.push_back(A);
	A = U[1]/(k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta));
	V1.push_back(A);
	A = -U[4];
	V1.push_back(A);
	A = -U[3];
	V1.push_back(A);
	A = (U[0] + U[2])/(k0*q0 - sqrt(k0*k0 - Mpiz*Mpiz)*sqrt(q0*q0 + Q2)*cos(theta));
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
	V.push_back(Value);
	Value = -178.45/((1.0+Q2)*(1.0+0.3457*Q2*Q2-0.087*Q2*Q2*Q2+0.00806*Q2*Q2*Q2*Q2));
	V.push_back(Value);
	Value = -339.06/((1.0+Q2)*(1.0+0.3481*Q2*Q2-0.0854*Q2*Q2*Q2+0.00758*Q2*Q2*Q2*Q2));
	V.push_back(Value);
	
	return V;
}

int main(int argc, char **argv) // W, Q2, cos(theta) generation
{
	srand(time(NULL));
	
	vector<complex<double>> Amp1, Amp3;
	vector<double> Amp2;
	
	double W, Q2, theta, k; k = 1000*Mpip; 
	
	W = fRand(1, 2); 
	Q2 = fRand(0, 4); 
	theta = fRand(0, M_PI); 
	
	Amp1 = Amplitudes_Fn(W, Q2, theta);
	Amp3 = Amplitudes_Fp(W, Q2, theta);
	Amp2 = Amplitudes_A(Q2);
	
	cout << "-----------------------------------------------" << endl;
	cout << "\tKinematic values\n\nW\t\t=\t" << W << "\tGeV\nQ2\t\t=\t" << Q2 << "\tGeV^2\nCos(Theta)\t=\t" << cos(theta) << endl;
	cout << "-----------------------------------------------" << endl;
	cout << "\tAmplitudes\n" << endl;
	cout << "(1) CGLN Amplitudes for Pi^{+}n final system:\n" << endl;
	cout << "\tF1\t=\t" << k*Amp1[0].real() << "\t+\t" << k*Amp1[0].imag() << endl;	
	cout << "\tF2\t=\t" << k*Amp1[1].real() << "\t+\t" << k*Amp1[1].imag() << endl;
	cout << "\tF3\t=\t" << k*Amp1[2].real() << "\t+\t" << k*Amp1[2].imag() << endl;
	cout << "\tF4\t=\t" << k*Amp1[3].real() << "\t+\t" << k*Amp1[3].imag() << endl;
	cout << "\tF5\t=\t" << k*Amp1[4].real() << "\t+\t" << k*Amp1[4].imag() << endl;
	cout << "\tF6\t=\t" << k*Amp1[5].real() << "\t+\t" << k*Amp1[5].imag() << endl;
	cout << "\n(2) CGLN Amplitudes for Pi^{0}p final system:\n" << endl;
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
	
	Amp1.clear(); Amp2.clear(); Amp3.clear();

	return 0;
}
