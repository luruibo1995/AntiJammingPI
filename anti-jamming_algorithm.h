#pragma once
#ifndef ANTI_JAMMING_ALGORITHM
#define ANTI_JAMMING_ALGORITHM
#include<fstream>
using namespace std;
const int Nmax = 16384;         //ÿ��ѭ������ĵ���
const double PI = 3.1416;
const int N = 4;               //����Ԫ
/***************�������㷨����*****************/
void PI_INV(ifstream& fin, ofstream& fout);
void frequency_filter(ifstream& fin, ofstream& fout);
void time_filter(ifstream& fin, ofstream& fout);
/***************�������㷨����*****************/

/************�����źŴ���������***************/
void FFT(double xreal[], double ximag[], int n);
void  IFFT(double xreal[], double ximag[], int n);
int Complex_Inv(int n, double ar[], double ai[]);//��������
inline void calculate_element(const double x1_r[], const double x1_i[], const double x2_r[], const double x2_i[], double& r, double& i);
inline void calculate_element(const double xr[], const double xi[], double& r, double& i);
/************�����źŴ���������***************/
#endif
