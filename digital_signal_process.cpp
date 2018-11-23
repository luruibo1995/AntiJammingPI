#include"anti-jamming_algorithm.h"
#include<math.h>



inline void swap(double &a, double &b)
{
	double t;
	t = a;
	a = b;
	b = t;
}

void bitrp(double xreal[], double ximag[], int n)
{
	// 位反转置换 Bit-reversal Permutation
	int i, j, a, b, p;

	for (i = 1, p = 0; i < n; i *= 2)
	{
		p++;
	}
	for (i = 0; i < n; i++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if (b > i)
		{
			swap(xreal[i], xreal[b]);
			swap(ximag[i], ximag[b]);
		}
	}
}

void FFT(double xreal[], double ximag[], int n)
{
	// 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
	double wreal[Nmax / 2], wimag[Nmax / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, index1, index2;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = -2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
				timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
				ureal = xreal[index1];
				uimag = ximag[index1];
				xreal[index1] = ureal + treal;
				ximag[index1] = uimag + timag;
				xreal[index2] = ureal - treal;
				ximag[index2] = uimag - timag;
			}
		}
	}
}

void  IFFT(double xreal[], double ximag[], int n)
{
	// 快速傅立叶逆变换
	double wreal[Nmax / 2], wimag[Nmax / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, index1, index2;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = 2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
				timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
				ureal = xreal[index1];
				uimag = ximag[index1];
				xreal[index1] = ureal + treal;
				ximag[index1] = uimag + timag;
				xreal[index2] = ureal - treal;
				ximag[index2] = uimag - timag;
			}
		}
	}

	for (j = 0; j < n; j++)
	{
		xreal[j] /= n;
		ximag[j] /= n;
	}
}

int Complex_Inv(int n, double ar[], double ai[])
{
	int *is, *js, i, j, k, l, u, v, w;
	double p, q, s, t, d, b;
	is = new int[n];
	js = new int[n];
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				u = i * n + j;
				p = ar[u] * ar[u] + ai[u] * ai[u];
				if (p>d) { d = p; is[k] = i; js[k] = j; }
			}
		if (is[k] != k)
		{
			j = 0;
			while (j <= n - 1)
			{
				u = k * n + j; v = is[k] * n + j;
				t = ar[u]; ar[u] = ar[v]; ar[v] = t;
				t = ai[u]; ai[u] = ai[v]; ai[v] = t;
				j++;
			}
		}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k; v = i * n + js[k];
				t = ar[u]; ar[u] = ar[v]; ar[v] = t;
				t = ai[u]; ai[u] = ai[v]; ai[v] = t;
			}
		l = k * n + k;
		ar[l] = ar[l] / d; ai[l] = -ai[l] / d;
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j;
				p = ar[u] * ar[l]; q = ai[u] * ai[l];
				s = (ar[u] + ai[u])*(ar[l] + ai[l]);
				ar[u] = p - q; ai[u] = s - p - q;
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				v = i * n + k;
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = k * n + j;  w = i * n + j;
						p = ar[u] * ar[v]; q = ai[u] * ai[v];
						s = (ar[u] + ai[u])*(ar[v] + ai[v]);
						t = p - q; b = s - p - q;
						ar[w] = ar[w] - t;
						ai[w] = ai[w] - b;
					}
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k;
				p = ar[u] * ar[l]; q = ai[u] * ai[l];
				s = (ar[u] + ai[u])*(ar[l] + ai[l]);
				ar[u] = q - p; ai[u] = p + q - s;
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = js[k] * n + j;
				t = ar[u]; ar[u] = ar[v]; ar[v] = t;
				t = ai[u]; ai[u] = ai[v]; ai[v] = t;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k; v = i * n + is[k];
				t = ar[u]; ar[u] = ar[v]; ar[v] = t;
				t = ai[u]; ai[u] = ai[v]; ai[v] = t;
			}
	}
	delete[] js; delete[] is;
	return(1);
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