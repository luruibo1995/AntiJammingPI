#pragma once
#ifndef ANTI_JAMMING_ALGORITHM
#define ANTI_JAMMING_ALGORITHM
#include<fstream>
using namespace std;
const int Nmax = 16384;         //每次循环计算的点数
const double PI = 3.1416;
const int N = 4;               //四阵元
/***************抗干扰算法声明*****************/
void PI_INV(ifstream& fin, ofstream& fout);
void frequency_filter(ifstream& fin, ofstream& fout);
void time_filter(ifstream& fin, ofstream& fout);
/***************抗干扰算法声明*****************/

/************数字信号处理函数声明***************/
void FFT(double xreal[], double ximag[], int n);
void  IFFT(double xreal[], double ximag[], int n);
int Complex_Inv(int n, double ar[], double ai[]);//仿真求逆
inline void calculate_element(const double x1_r[], const double x1_i[], const double x2_r[], const double x2_i[], double& r, double& i);
inline void calculate_element(const double xr[], const double xi[], double& r, double& i);
/************数字信号处理函数声明***************/
#endif
