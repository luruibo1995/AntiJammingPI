#include"anti-jamming_algorithm.h"
#include<iostream>
#include<fstream>
#include<time.h>
using namespace std;



void main()
{
	const int fs = 62000000;//������
	int loop;  //ѭ������
	int t = 40;//����ʱ�䣬��λΪ�루s��
	loop = fs * t / Nmax;

	ifstream fread("G:\\CW.bin", ios::in | ios::binary);
	ofstream fwrite("G:\\yout.dat", ios::out | ios::binary);

	/*clock_t start, finish;
	double totaltime;
	start = clock();
*/

	for (int i = 0; i < 1; i++)
	{
		//frequency_filter(fread, fwrite);
		PI_INV(fread,fwrite); 
		//PI_INV_2(fread, fwrite);
		time_filter(fread,fwrite);
	}

	

	//finish = clock();
	//totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "\n�˳��������ʱ��Ϊ" << totaltime << "�룡" << endl;

	fread.close();
    fwrite.close();


	system("pause");
}
