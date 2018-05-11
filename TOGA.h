# TOGA
A famous lightning location method
#pragma once
#include <vector>
#include <string>
#include <complex>
#include <math.h>
#include <time.h>
#include <iostream>

using namespace std;

struct Shuchu {
	unsigned int year;
	unsigned char month;
	unsigned char day;
	unsigned char hour;
	unsigned char minute;
	unsigned char second;
	unsigned long microsecond;
	short int ID;
	double Peaks;
	double W_mod;
	double Eng;
	double Vg;
	int flag;
	complex<double> *fft=new complex<double> [257]();
	double *phase = new double[256]();
	double *wave=new double [1020]();
};

class TOGA
{
public:
	TOGA();
	~TOGA();
	struct Shuchu Read(unsigned char *buf);
	complex<double>* wavedect(double *Wave, complex<double> *Complex, int Id_max);
	vector<vector<int> >find_Peaks(double *Wave, int len);
	int find_next(int point, vector<int> Ind, double *Wave);
	complex<double>* fft(complex<double> *Wave, int Nf, int flag);
	double trapz(double *Wave, int start, int end);
	int zero_crossing(int* time, double*  Wave);
	int lagrange(int *Y, double *X, double x,int N);
	struct Shuchu calculate(struct Shuchu Out);
	double* linspace(double a, double b, double n);
	double* cal_phase(complex<double>* FFT);
	double rad2deg(double point);
	double deg2rad(double point);
	double angle(complex<double> point);
	double* lagrpol(vector<double> X, vector<double> Y, double *x);
	double TOGA::Optimize(double *X, double *Y);
private:
	FILE *file;
};
