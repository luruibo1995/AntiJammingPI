#include"pch.h"
#include"anti-jamming_algorithm.h"
#include<fstream>
#include<iostream>
#include<time.h>
using namespace std;

inline short pow_2(int n);

char* PI_INV(char* buf)
{
	/*****************************************/
	clock_t start, finish;
	double totaltime;
	/*****************************************/

	//short x[Nmax * 4];
	double* x1_r = new double[Nmax];
	double* x1_i = new double[Nmax]();
	double* x2_r = new double[Nmax];
	double* x2_i = new double[Nmax]();
	double* x3_r = new double[Nmax];
	double* x3_i = new double[Nmax]();
	double* x4_r = new double[Nmax];
	double* x4_i = new double[Nmax]();

	short* x = (short*)buf;

	//start = clock();	
	for (int i = 0; i < (Nmax * 4); i++)
	{
		int n = (i / 4);
		if (i % 4 == 0)
			x1_r[n] = (double)x[i];
		else if (i % 4 == 1)
			x2_r[n] = (double)x[i];
		else if (i % 4 == 2)
			x3_r[n] = (double)x[i];
		else x4_r[n] = (double)x[i];
	}

	for (int i = 0; i < 10; i++)
	{
		cout << x1_r[i] << endl;
	}
	/*finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;
*/
//start = clock();
/*********************************************************希尔伯特变换*********************************************************/
	FFT(x1_r, x1_i, Nmax);
	FFT(x2_r, x2_i, Nmax);
	FFT(x3_r, x3_i, Nmax);
	FFT(x4_r, x4_i, Nmax);

	for (int i = 0; i < Nmax / 2; i++)
	{
		x1_r[i] = 2 * x1_r[i];
		x1_i[i] = 2 * x1_i[i];
		x2_r[i] = 2 * x2_r[i];
		x2_i[i] = 2 * x2_i[i];
		x3_r[i] = 2 * x3_r[i];
		x3_i[i] = 2 * x3_i[i];
		x4_r[i] = 2 * x4_r[i];
		x4_i[i] = 2 * x4_i[i];
		x1_r[i + 8192] = 0;
		x1_i[i + 8192] = 0;
		x2_r[i + 8192] = 0;
		x2_i[i + 8192] = 0;
		x3_r[i + 8192] = 0;
		x3_i[i + 8192] = 0;
		x4_r[i + 8192] = 0;
		x4_i[i + 8192] = 0;
	}
	IFFT(x1_r, x1_i, Nmax);
	IFFT(x2_r, x2_i, Nmax);
	IFFT(x3_r, x3_i, Nmax);
	IFFT(x4_r, x4_i, Nmax);
	/*********************************************************希尔伯特变换*********************************************************/
	/*finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;
*/
//start = clock();
/*********************************************************求协方差矩阵*********************************************************/
	double Rr[4][4], Ri[4][4];
	double r, i;

	calculate_element(x1_r, x1_i, r, i);
	Rr[0][0] = r; Ri[0][0] = i;
	calculate_element(x1_r, x1_i, x2_r, x2_i, r, i);
	Rr[0][1] = r; Ri[0][1] = i; Rr[1][0] = r; Ri[1][0] = -1 * i;
	calculate_element(x1_r, x1_i, x3_r, x3_i, r, i);
	Rr[0][2] = r; Ri[0][2] = i; Rr[2][0] = r; Ri[2][0] = -1 * i;
	calculate_element(x1_r, x1_i, x4_r, x4_i, r, i);
	Rr[0][3] = r; Ri[0][3] = i; Rr[3][0] = r; Ri[3][0] = -1 * i;
	calculate_element(x2_r, x2_i, r, i);
	Rr[1][1] = r; Ri[1][1] = i;
	calculate_element(x2_r, x2_i, x3_r, x3_i, r, i);
	Rr[1][2] = r; Ri[1][2] = i; Rr[2][1] = r; Ri[2][1] = -1 * i;
	calculate_element(x2_r, x2_i, x4_r, x4_i, r, i);
	Rr[1][3] = r; Ri[1][3] = i; Rr[3][1] = r; Ri[3][1] = -1 * i;
	calculate_element(x3_r, x3_i, r, i);
	Rr[2][2] = r; Ri[2][2] = i;
	calculate_element(x3_r, x3_i, x4_r, x4_i, r, i);
	Rr[2][3] = r; Ri[2][3] = i; Rr[3][2] = r; Ri[3][2] = -1 * i;
	calculate_element(x4_r, x4_i, r, i);
	Rr[3][3] = r; Ri[3][3] = i;
	/*********************************************************求协方差矩阵*********************************************************/
	//finish = clock();
	//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;


	/*********************************************************求矩阵的逆*********************************************************/
	int k = Complex_Inv(N, &Rr[0][0], &Ri[0][0]);
	/*********************************************************求矩阵的逆*********************************************************/

	/*********************************************************求权值wopt*********************************************************/
	double Wopt_r[N], Wopt_i[N];
	for (int i = 0; i < N; i++)
	{
		Wopt_r[i] = (Rr[i][0] * Rr[0][0] + Ri[i][0] * Ri[0][0]) / (Rr[0][0] * Rr[0][0] + Ri[0][0] * Ri[0][0]);
		Wopt_i[i] = (Ri[i][0] * Rr[0][0] - Rr[i][0] * Ri[0][0]) / (Rr[0][0] * Rr[0][0] + Ri[0][0] * Ri[0][0]);
	}
	/*********************************************************求权值wopt*********************************************************/

	/*********************************************************求输出yout*********************************************************/
	double yout[Nmax];//这里取输出的实部，所以这里只有输出信号的实部
	for (int i = 0; i < Nmax; i++)
	{
		yout[i] = (Wopt_r[0] * x1_r[i]) + (Wopt_i[0] * x1_i[i]) + (Wopt_r[1] * x2_r[i]) + (Wopt_i[1] * x2_i[i]) + \
			(Wopt_r[2] * x3_r[i]) + (Wopt_i[2] * x3_i[i]) + (Wopt_r[3] * x4_r[i]) + (Wopt_i[3] * x4_i[i]);
	}
	/*********************************************************求输出yout*********************************************************/
	/*finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n此程序的运行时间为" << totaltime << "秒！" << endl;
*/
/*start = clock();*/

/*******************************************************输出.dat文件*********************************************************/
//unsigned data = 0;
//for (int i = 0; i < (Nmax / 8); i++)
//{
//	data = 0;
//	for (int j = 0; j < 8; j++)
//	{
//		if (yout[i * 8 + j] < 0)
//			data = data + pow_2(7 - j);
//	}
//	fout.write((char*)&data, 1);
//}
	unsigned data[Nmax / 16] = {0};
	for (int i = 0; i < (Nmax / 16); i++)
	{
		for (int j = 0; j < 16; j++)
		{
			if (yout[i * 16 + j] < 0)
				data[i] = data[i] + pow_2(15 - j);
		}
	}
	/*******************************************************输出.dat文件*********************************************************/
	
	/**************************测试***************************/
	//ofstream fout1("G:\\ch1.dat", ios::out | ios::binary);
	//ofstream fout2("G:\\ch2.dat", ios::out | ios::binary);
	//ofstream fout3("G:\\ch3.dat", ios::out | ios::binary);
	//ofstream fout4("G:\\ch4.dat", ios::out | ios::binary);
	//ofstream fout5("G:\\ch5.dat", ios::out | ios::binary);
	//ofstream fout6("G:\\ch6.dat", ios::out | ios::binary);
	//ofstream fout7("G:\\ch7.dat", ios::out | ios::binary);
	//ofstream fout8("G:\\ch8.dat", ios::out | ios::binary);
	//ofstream fout9("G:\\yout1.dat", ios::out | ios::binary);

	//for (int i = 0; i < Nmax; i++)
	//{
	//	//fout1 << x1[i];
	//	//fout2 << x2[i];
	//	//fout3 << x3[i];
	//	//fout4 << x4[i];
	//	fout1.write((const char*)&x1_r[i], 8);
	//	fout2.write((const char*)&x1_i[i], 8);
	//	fout3.write((const char*)&x2_r[i], 8);
	//	fout4.write((const char*)&x2_i[i], 8);
	//	fout5.write((const char*)&x3_r[i], 8);
	//	fout6.write((const char*)&x3_i[i], 8);
	//	fout7.write((const char*)&x4_r[i], 8);
	//	fout8.write((const char*)&x4_i[i], 8);
	//	fout9.write((const char*)&yout[i], 8);
	//}
	//fout1.close();
	//fout2.close();
	//fout3.close();
	//fout4.close();
	//fout5.close();
	//fout6.close();
	//fout7.close();
	//fout8.close();
	//fout9.close();
	/**************************测试***************************/
	delete[] x1_r;
	delete[] x1_i;
	delete[] x2_r;
	delete[] x2_i;
	delete[] x3_r;
	delete[] x3_i;
	delete[] x4_r;
	delete[] x4_i;

	return (char*)data;
}


inline void calculate_element(const double x1_r[], const double x1_i[], const double x2_r[], const double x2_i[], double& r, double& i)
{
	r = 0;
	i = 0;
	for (int j = 0; j < Nmax; j++)
	{
		r = r + x1_r[j] * x2_r[j] + x1_i[j] * x2_i[j];
		i = i + x2_r[j] * x1_i[j] - x1_r[j] * x2_i[j];
	}
}
inline void calculate_element(const double xr[], const double xi[], double& r, double& i)
{
	r = 0;
	i = 0;
	for (int j = 0; j < Nmax; j++)
	{
		r = r + xr[j] * xr[j] + xi[j] * xi[j];
	}
}

inline short pow_2(int n)
{
	short x = 1;
	for (int i = 0; i < n; i++)
	{
		x = x << 1;
	}
	return x;
}