#include<iostream>
#include<fstream>

using namespace std;

#define Nmax 16384*2

//void calculate_element(const double x1_r[], const double x1_i[], const double x2_r[], const double x2_i[], double& r, double& i);
//void calculate_element(const double xr[], const double xi[], double& r, double& i);

void calculate_element(double& y,const double r[],const double i[], int n);
void what();


int main()
{
	short x = 3;
	char c;
	
	c = (char)x;
	
	

	

	cout << x << c << endl;
	what();
	system("pause");
	return 0;
}

void calculate_element(double& y, const double r[], const double i[],int n)
{
	for(int j = 0;j < n;j++)
	y = r[j] + i[j];
}

//void calculate_element(double x1_r[], double x1_i[], double x2_r[], double x2_i[], double& r, double& i)
//{
//	r = 0;
//	i = 0;
//	for (int j = 0; j < Nmax; j++)
//	{
//		r = r + x1_r[j] * x2_r[j] + x1_i[j] * x2_i[j];
//		i = i + x2_r[j] * x1_i[j] - x1_r[j] * x2_i[j];
//	}
//}
//void calculate_element(double xr[], double xi[], double& r, double& i)
//{
//	r = 0;
//	i = 0;
//	for (int j = 0; j < Nmax; j++)
//	{
//		r = r + xr[j] * xr[j] + xi[j] * xi[j];
//	}
//}

void what()
{
	using namespace std;
	cout << "zjl sb" << endl;
}