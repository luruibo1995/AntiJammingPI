#include"anti-jamming_algorithm.h"
#include<fstream>
using namespace std;

inline short pow_2(int n);

void frequency_filter(ifstream& fin, ofstream& fout)
{


	short x[Nmax * 4];
	double* x1_r = new double[Nmax];
	double* x1_i = new double[Nmax]();
	double* h1_r = new double[Nmax]();
	double* h1_i = new double[Nmax]();
	double  fy_r[Nmax] = { 0 }, fy_i[Nmax] = { 0 };

	for (int i = 0; i < (Nmax * 4); i++)
	{
		fin.read((char*)&x[i], 2);
	}

	for (int i = 0; i < (Nmax * 4); i++)
	{
		int n = (i / 4);
		/*cout << i << "		" << n << endl;*/
		if (i % 4 == 0)
			x1_r[n] = (double)x[i];
	}
	FFT(x1_r, x1_i, Nmax);
	for (int i = 0; i < Nmax / 2; i++)
	{
		h1_r[i] = 2 * x1_r[i];
		h1_i[i] = 2 * x1_i[i];
	}
	for (int i = 0; i < Nmax / 2; i++)
	{
		double n = h1_r[i] * h1_r[i] + h1_i[i] * h1_i[i];
		if (n < 12841195761)
		{
			fy_r[i] = h1_r[i];
			fy_i[i] = h1_i[i];
		}
	}
	IFFT(fy_r, fy_i, Nmax);
	/*********************************************写文件***************************************************/
	unsigned data = 0;
	for (int i = 0; i < (Nmax / 8); i++)
	{
		data = 0;
		for (int j = 0; j < 8; j++)
		{
			if (fy_r[i * 8 + j] < 0)
				data = data + pow_2(7 - j);
		}
		fout.write((char*)&data, 1);
	}
	/*********************************************写文件***************************************************/

	delete[] x1_r;
	delete[] x1_i;
	delete[] h1_r;
	delete[] h1_i;

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