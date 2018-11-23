#include<iostream>
#include"anti-jamming_algorithm.h"


inline short pow_2(int n)
{
	short x = 1;
	for (int i = 0; i < n; i++)
	{
		x = x << 1;
	}
	return x;
}

void t_f_filter(ifstream& fin, ofstream& fout)
{
	short x[Nmax * 4 + 256];
	double** xr = new double*[64]; double** xi = new double*[64];
	/*double** xr[64]; double** xi[64];*/

	for (int i = 0; i < 64; i++)
	{
		xr[i] = new double[Nmax];
		xi[i] = new double[Nmax]();
	}

	for (int i = 0; i < (Nmax * 4 + 256); i++)
	{
		fin.read((char*)&x[i], 2);
	}

	for (int j = 63; j >= 0; j--)
	{
		for (int i = 0; i < (Nmax * 4); i++)
		{
			int n = (i / 4);
			if (i % 4 == 0)
			{
				xr[j][n] = (double)x[i + (63 - j) * 4];
			}
		}
	}
	for (int i = 0; i < 64; i++)
		FFT(xr[i], xi[i], Nmax);

	for (int j = 0; j < 64; j++)
		for (int i = 0; i < Nmax / 2; i++)
		{
			xr[j][i] = 2 * xr[j][i];
			xi[j][i] = 2 * xi[j][i];
			xr[j][i + 8192] = 0;
			xi[j][i + 8192] = 0;
		}
	for(int j = 0;j < 64; j++)
		for (int i = 0; i < Nmax / 2; i++)
		{
			double n = xr[j][i] * xr[j][i] + xi[j][i] * xi[j][i];
			if (n > 12841195761)
			{
				xr[j][i] = 0;
				xi[j][i] = 0;
			}
		}

	for (int i = 0; i < 64; i++)
		IFFT(xr[i], xi[i], Nmax);

	double Rr[64][64], Ri[64][64] = { 0 };
	double r = 0, i = 0;
	for (int n = 1; n < 64; n++)
	{

		calculate_element(xr[0], xi[0], xr[n], xi[n], r, i); r = r / Nmax; i = i / Nmax;
		for (int m = 0; m < 64 - n; m++)
		{
			Rr[m][n + m] = r; Ri[m][n + m] = i;
			Rr[n + m][m] = r; Ri[n + m][m] = -1 * i;
		}
	}
	calculate_element(xr[0], xi[0], xr[0], xi[0], r, i); r = r / Nmax; i = i / Nmax;
	for (int j = 0; j < 64; j++)
	{
		Rr[j][j] = r;
		Ri[j][j] = i;
	}
	double rx_r[64], rx_i[64];
	for (int i = 0; i < 64; i++)
	{
		rx_r[i] = Rr[0][i];
		rx_i[i] = Ri[i][0];
	}
	int k = 0;
	k = Complex_Inv(64, &Rr[0][0], &Ri[0][0]);

	double Wopt_r[64], Wopt_i[64];
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			Wopt_r[i] = (Rr[i][j] * rx_r[j] + Ri[i][j] * rx_i[j]);
			Wopt_i[i] = (Ri[i][j] * rx_r[j] - Rr[i][j] * rx_i[j]);
		}
	}
	double yout[Nmax];
	for (int i = 0; i < Nmax; i++)
	{
		for (int j = 0; j < 64; j++)
		{
			yout[i] = Wopt_r[j] * xr[j][i] + Wopt_i[j] * xi[j][i];
		}
	}
	/*******************************************************输出.dat文件*********************************************************/
	unsigned data = 0;
	for (int i = 0; i < (Nmax / 8); i++)
	{
		data = 0;
		for (int j = 0; j < 8; j++)
		{
			if (yout[i * 8 + j] < 0)
				data = data + pow_2(7 - j);
		}
		fout.write((char*)&data, 1);
	}
	/*******************************************************输出.dat文件*********************************************************/
	for (int i = 0; i < 64; i++)
	{
		delete[] xr[i];
		delete[] xi[i];
	}
	delete[] xr;
	delete[] xi;

}


