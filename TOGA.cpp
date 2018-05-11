# TOGA
A famous lightning location method

#include <algorithm>
#include <cmath>
#include <iostream>
#include"TOGA.h"
using namespace std;
extern const double pi = 3.14159265358979323846;
extern const int SN = 2048;
extern const int Nf = 256;
extern const int Nf2 = 1000;
TOGA::TOGA()
{
	cout << "TOGA Constructors" << endl;
}

TOGA::~TOGA()
{
	if (file)
	{
		fclose(file);
	}
}



struct Shuchu TOGA::Read(unsigned char *buf)
{
	Shuchu Out;
	complex<double> *Complex_data=new complex<double>[(SN - 8) / 2]();
	double Max = 0;
	double Min = 0;
	unsigned int Id_max = 0;
	unsigned int Id_min = 0;
	int temp1;
	int temp2;
	const double uncorrect = -4.0;
	const double MUL = 256.0;
	Out.year = buf[4087] + 2000;
	Out.month = buf[4086];	
	Out.day = buf[4089];
	Out.hour = buf[4088];
	Out.minute = buf[4091];
	Out.second = buf[4090];
	unsigned long us = buf[4092] + (buf[4093] << 8) + (buf[4094] << 16);
	unsigned long Zhengshu = us / 10;
	char Namiao = us % 10;
	short int SID = int((buf[5] << 2) + (buf[4] >> 6));
	Out.ID = SID;
	for (int i = 1; i < 2040; i += 2)
	{
		temp1 = buf[2 * i + 6];
		temp2 = buf[2 * i + 7];
		Out.wave[i / 2] = ((temp1 + int(temp2) * 256) - 2047)*10.0 / 4095.0;
		if (Out.wave[i / 2] > Max)
		{
			Max = Out.wave[i / 2];
			Id_max = i / 2;
		}
		if (Out.wave[i / 2] < Min)
		{
			Min = Out.wave[i / 2];
			Id_min = i / 2;
			if (Min < uncorrect)
			{
				Out.wave[i / 2] = 0;
			}
		}
		Complex_data[i / 2] = complex<double>(Out.wave[i / 2], 0);		
	}
	Out.fft = wavedect(Out.wave, Complex_data, int(Id_max));
	delete Complex_data;
	Out.Eng = abs(trapz(Out.wave, 0, 1019));
	Out.flag = Out.fft[256].real();
	int delta = Id_max - 120;
	Out.microsecond = (Zhengshu + delta) * 10 + Namiao;
	return Out;
}

complex<double>* TOGA::wavedect(double *Wave, complex<double> *Complex, int Id_max) 
{
	int Flag;
	const int flag = 1;
	const int Len = (SN - 8) / 2;
	vector<vector <int> > Ind;
	Ind = find_Peaks(Wave, Len);
	int next_min_point = find_next(Id_max, Ind[1], Wave); 
	int next_max_point = find_next(next_min_point, Ind[0], Wave);
	complex<double>* FFT = fft(Complex, Nf, flag);
	double *amp = new double[Nf]();
	for (int i = 0; i < Nf; i++)
	{		
		amp[i] = abs(FFT[i]);		
	}
	double Inte1 = trapz(amp, 1, 6);
	double Inte2 = trapz(amp, 7, 16);
	double Inte3 = trapz(amp, 17, 25);
	if (next_min_point == Id_max)
	{
		Flag = 2;
	}
	else if (next_min_point == next_max_point)
	{
		Flag = 2;
	}
	else if (next_min_point != Id_max)
	{
		int time1[2] = { Id_max, next_min_point };
		int Zero_time1 = zero_crossing(time1, Wave);			
		int time2[2] = { next_min_point, next_max_point };
		int Zero_time2 = zero_crossing(time2, Wave);
		if ((Wave[Id_max] > 4.5) || (Inte1 / (Inte1 + Inte2 + Inte3) < 0.4))///近处波
		{
			Flag = 1;
		}
		else if ((Inte1 / (Inte1 + Inte2 + Inte3) >= 0.4))///远处波
		{
			Flag = 2;
		}
	}
	Ind.clear();
	delete amp;
	FFT[Nf] = Flag;
	return FFT;
}
vector<vector<int> > TOGA::find_Peaks(double *Wave, int len)
{
	double *p = Wave;//引入矩阵
	double const Max = 2;
	double const Min = -2;
	double *Diffv = new double[len - 1]();
	double *Trend = new double[len - 1]();
	for (int i = 1; i < len; i++)
	{
		Diffv[i - 1] = p[i] - p[i - 1];

		if (Diffv[i - 1] > 0)
		{
			Trend[i - 1] = 1;
		}
		else if (Diffv[i - 1] == 0)
		{
			Trend[i - 1] = 0;
		}
		else if (Diffv[i - 1] < 0)
		{
			Trend[i - 1] = -1;
		}
	}
	for (int i = len - 2; i >= 0; i--)
	{

		if ((Trend[i] == 0) && (i == len - 2))
		{
			Trend[i] = 1;
		}
		else if (Trend[i] == 0)
		{
			if (Trend[i + 1] >= 0)
			{
				Trend[i] = 1;
			}
			else
			{
				Trend[i] = -1;
			}
		}
	}
	double *R = new double[len - 2];
	vector<vector<int> >Ind(2);//二维数组
	for (int i = 1; i < len - 1; i++)
	{
		R[i - 1] = Trend[i] - Trend[i - 1];
		if (R[i - 1] == Max)
		{
			Ind[1].push_back(i);
		}
		if (R[i - 1] == Min)
		{
			Ind[0].push_back(i);
		}
	}
	delete R;
	delete Diffv;
	delete Trend;
	return Ind;
}

int TOGA::find_next(int point, vector<int> Ind, double *Wave)
{
	vector<int> temp;
	for (int i = 0; i < Ind.size(); i++)
	{
		if (Ind[i] > point)
		{
			if ((Wave[point] > 0) && (Wave[Ind[i]] < 0))//判断波谷
			{
				temp.push_back(Ind[i] - point);
			}
			else if ((Wave[point] < 0) && (Wave[Ind[i]] > 0))//判断波峰
			{
				temp.push_back(Ind[i] - point);
			}
		}
	}
	if (temp.empty())
	{
		return point;
	}
	else
	{
		return point + temp[0];
	}
	temp.clear();
}


complex<double>* TOGA::fft(complex<double> *Wave, int Nf, int flag)
{
	const double a = 2.0;
	complex<double>  *temp = new complex<double>[Nf+1];
	for (int i = 0; i < Nf; i++)
	{
		temp[i] = 0;
		for (int j = 0; j < Nf; j++)
		{
			complex<double> wn = complex<double>(cos(a*pi / Nf*i*j), -sin(flag*a*pi / Nf*i*j));
			temp[i] += Wave[j] * wn;			
		}	
	}
	if (flag = 1)
	{
		return temp;
	}		
}

double TOGA::trapz(double *Wave, int start, int end)
{
	const int n = 1;//dx=1;
	const int h = (end - start+1) / n;///分成h部分
	double w = 0.0;
	for (int i = 0; i <h - 3; i++)
	{
		w += Wave[start + (i + 1) ];
		
	}
	return (Wave[start] + Wave[end-1] + 2 * w) / 2;
}

int TOGA::zero_crossing(int *time, double *Wave)
{
	double *amp = new double[time[1] - time[0]+ 1]();
	const int N = time[1] - time[0] + 1;//间隔
	for (int t = time[0]; t <=time[1]; t++)
	{
		amp[t - time[0]] = Wave[t];	
	}
	int zero_time = lagrange(time, amp, 0.0, N);
	delete amp;
	return zero_time;
}

int TOGA::lagrange(int *Y, double *X, double x,int N)
{
	for (int i = 0; i < N - 1; i++)
	{
		int temp = (Y[0] + i);
		if ((X[i] < 0.0) && (X[i + 1]>0.0))
		{
			return temp + 1;
			break;
		}
		if ((X[i] > 0.0) && (X[i + 1] < 0.0))
		{
			return temp;
			break;
		}
	}
}

struct Shuchu TOGA::calculate(struct Shuchu Out)
{
	if (Out.flag == 2)
	{
		double* f = linspace(0, 1000000, Nf);
		vector<double> f_sample;//取f中4500到24000的
		vector<double> phase_sample; //取f中4500到24000对应的相位
		double* phase = cal_phase(Out.fft);
		int n = 0;
		for (int i = 0; i < Nf; i++)
		{				
			if ((f[i] >= 3500) && (f[i] < 24000))
			{
				f_sample.push_back(f[i]);
				phase_sample.push_back(phase[i]);
				n++;
			}
		}
		double *f_sample1 = new double[n]();
		double *phase_sample1 = new double[n]();
		for (int i = 0; i < n; i++)
		{
			f_sample1[i] = f_sample[i];
			phase_sample1[i] = phase_sample[i];
		}
		//double *f_match = linspace(6000, 18000, Nf2);
		//double *y = lagrpol(f_sample, phase_sample, f_match);
		Out.W_mod = Optimize(f_sample1, phase_sample1);
		double W0 = 1670;
		const double c = 2.997925*1e8;
		Out.Vg = c*(sqrt(1 - (pow(W0, 2) / pow(Out.W_mod, 2))));
		//delete f_match;
		delete f;
		f_sample.clear();
		phase_sample.clear();
		//delete y;
		delete phase;
		delete f_sample1;
		delete phase_sample1;
	}
	if (Out.flag == 1)
	{
		Out.Vg= 2.997925*1e8;
	}
	return Out;
}

double* TOGA::lagrpol(vector<double> X, vector<double> Y, double *x)
{
	double *Out = new double[Nf2]();
	for (int i = 0; i < Nf2; i++)
	{
		float y = 0.0;
		float temp = x[i];
		for (int k = 0; k < X.size(); k++)
		{			
			float t = 1.0;
			for (int j = 0; j < X.size(); j++)
			{
				if (k != j)
				{
					t = t*((temp - X[j]) / (X[k] - X[j]));
				}
			}
			y = y + Y[k] * t;
		}
		Out[i] = y;			
	}
	return Out;
}

double* TOGA::linspace(double a, double b, double n)
{
	double* array=new double[n]();  //申请一个空间
	double step = (b - a) / (n - 1);
	for (int i=0;i<n;i++)
	{
		array[i]=a;//a存到array中		
		a += step;			
	}
	return array;
}

double* TOGA::cal_phase(complex<double>*  FFT)
{
	double *point = new double[Nf]();
	double *deg = new double[Nf]();
	double *phase = new double[Nf]();
	double *temp = new double[Nf]();
	for (int i = 0; i < Nf; i++)
	{				
		point[i] = angle(FFT[i]);
		deg[i] = rad2deg(point[i] - pi / 2);
		if (deg[i] <= -180.)
		{
			deg[i] = 360. + deg[i];
		}
		else if ((deg[i] <= 0.0) && (deg[i]>-180.0))
		{
			deg[i] = abs(deg[i]);
		}
		temp[i] = deg2rad(deg[i]);
		if (deg[i] == 180.0)
		{
			temp[i] = pi;
		}
		phase[i] = -acos(temp[i] / pi) / pi;		
	}
	delete temp;
	delete point;
	delete deg;
	return phase;
}

double TOGA::rad2deg(double point)
{
	return point*180. / pi;
}

double TOGA::deg2rad(double point)
{
	return point*pi / 180.;
}

double TOGA::angle(complex<double> point)
{
	double Angle = 0;
	if (point.real() >= 0 && point.imag() >= 0)
	{
		Angle = atan(point.imag() / point.real());
	}
	if (point.imag() > 0 && point.real() < 0)
	{
		Angle = atan(point.imag() / point.real()) + pi;
	}
	if (point.real() < 0 && point.imag() < 0)
	{
		Angle = atan(point.imag() / point.real()) - pi;
	}
	if (point.imag() < 0 && point.real() > 0)
	{
		Angle = atan(point.imag() / point.real());
	}
	return Angle;
}

double TOGA::Optimize(double *X, double *Y)
{
	
	vector<vector<int> > peaks = find_Peaks(Y, Nf2);
	return X[peaks[1][0]];
}
