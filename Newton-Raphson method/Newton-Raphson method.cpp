// Newton-Raphson method.cpp: ���������� ����� ����� ��� ����������� ����������.
//

#include "stdafx.h"
#include "math.h"
#include <iostream>
using namespace std;

double f(double x, double y, double s, double dk0, double dk1) //���������� �������� ������� �������
{
	double X, Y;
	X = x + s * dk0;
	Y = y + s * dk1;
	return X * X - X * Y + Y * Y + 9 * X - 6 * Y + 20;
}

double dfx(double x, double y, double s, double dk0, double dk1) //���������� �������� ������ ����������� �� x
{
	double dx = 0.001;
	return (f(x + dx, y, s, dk0, dk1) - f(x, y, s, dk0, dk1)) / dx;
}

double dfy(double x, double y, double s, double dk0, double dk1) //���������� �������� ������ ����������� �� y
{
	double dy = 0.001;
	return (f(x, y + dy, s, dk0, dk1) - f(x, y, s, dk0, dk1)) / dy;
}

double d2fx(double x, double y, double s, double dk0, double dk1) // ���������� �������� ������ ����������� �� x
{
	double dx = 0.001;
	return (f(x + dx, y, s, dk0, dk1) - 2 * f(x, y, s, dk0, dk1) + f(x - dx, y, s, dk0, dk1)) / (dx*dx);
}

double d2fy(double x, double y, double s, double dk0, double dk1) // ���������� �������� ������ ����������� �� y
{
	double dy = 0.001;
	return (f(x, y + dy, s, dk0, dk1) - 2 * f(x, y, s, dk0, dk1) + f(x, y - dy, s, dk0, dk1)) / (dy*dy);
}

double d2fxy(double x, double y, double s, double dk0, double dk1) // ���������� �������� ������ ����������� �� xy
{
	double dx = 0.001;
	double dy = 0.001;
	return (dfy(x + dx, y, s, dk0, dk1) - dfy(x, y, s, dk0, dk1)) / dx;
}

double dfs(double x, double y, double s, double dk0, double dk1) // ���������� �������� ����������� �� s
{
	double ds = 0.001;
	return (f(x, y, s + ds, dk0, dk1) - f(x, y, s, dk0, dk1)) / ds;
}

double d2fs(double x, double y, double s, double dk0, double dk1)  // ���������� �������� ������ ����������� �� s
{
	double ds = 0.001;
	return (f(x, y, s + ds, dk0, dk1) - 2 * f(x, y, s, dk0, dk1) + f(x, y, s - ds, dk0, dk1)) / (ds*ds);
}

double d3fs(double x, double y, double s, double dk0, double dk1)  // ���������� �������� ������ ����������� �� s
{
	double ds = 0.001;
	return (dfs(x, y, s + ds, dk0, dk1) - 2 * dfs(x, y, s, dk0, dk1) + dfs(x, y, s - ds, dk0, dk1)) / (ds*ds*ds);
}


int main()
{
	setlocale(0, "Rus");
	double x0, y0, x, y, k = 0;
	double a, b, s0, s, flag = 0;
	double eps1, eps2, dk[2], ff[2];
	double gesse[2][2], det_gesse, obr_gesse[2][2];
	cout << "� � � � �   � � � � � � � - � � � � � � �" << endl;
	cout << "f(x,y)=x*x-x*y+y*y+9*x-6*y+20 --> min" << endl;
	cout << "\n������� ��������� ����������� x0 � y0" << endl;
	cout << "x0= "; cin >> x0;
	cout << "y0= "; cin >> y0;
	cout << "������� ��������" << endl;
	cout << "eps1= "; cin >> eps1;
	cout << "eps2= "; cin >> eps2;
	cout << "������� ������� ������� [a;b] ��� ������ ���� s ������� �������:" << endl;
	cout << "a= "; cin >> a;
	cout << "b= "; cin >> b;
	dk[0] = 0;
	dk[1] = 0;
	cout << endl;
	system("pause");
	system("cls");
	cout << "� � � � �   � � � � � � � - � � � � � � �" << endl;
	cout << "f(x,y)=x*x-x*y+y*y+9*x-6*y+20 --> min" << endl;
	cout << "_________________________________________" << endl;
	cout << "�������� ������� � ����������� � ����� (x0,y0): " << endl;
	cout << "f(x,y) = " << f(x0, y0, 0, 0, 0) << endl;
	cout << "df(x,y)/dx = " << dfx(x0, y0, 0, 0, 0) << endl;
	cout << "df(x,y)/dy = " << dfy(x0, y0, 0, 0, 0) << endl;
	cout << "d^2f(x,y)/dxdy = " << d2fxy(x0, y0, 0, 0, 0) << endl;
	cout << "d^2f(x,y)/dx^2 = " << d2fx(x0, y0, 0, 0, 0) << endl;
	cout << "d^2f(x,y)/dy^2 = " << d2fy(x0, y0, 0, 0, 0) << endl;
	if (sqrt(dfx(x0, y0, 0, 0, 0)*dfx(x0, y0, 0, 0, 0) + dfy(x0, y0, 0, 0, 0)*dfy(x0, y0, 0, 0, 0)) <= eps1)
		cout << "������� � ����� x*=" << x0 << " y*=" << y0;
	else
		do
		{
			if (flag == 2)
			{
				x0 = x;
				y0 = y;
			}
			k++;
			cout << endl << "_____�������� �" << k << "_____";
			//������� �����
			gesse[0][0] = d2fx(x0, y0, 0, 0, 0);  gesse[0][1] = d2fxy(x0, y0, 0, 0, 0);
			gesse[1][0] = d2fxy(x0, y0, 0, 0, 0); gesse[1][1] = d2fy(x0, y0, 0, 0, 0);
			cout << endl << "***������� �����***" << endl;
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
					printf(" %7.3f", gesse[i][j]);
				cout << endl;
			}
			//���������� �������� ������� �����
			det_gesse = gesse[0][0] * gesse[1][1] - gesse[0][1] * gesse[1][0];
			cout << "������������ ������� ����� det=" << det_gesse << ">0" << endl;
			if (det_gesse == 0)
				printf("������� �����������, �������� ������� ����� ������!\n");
			else
			{
				obr_gesse[0][0] = gesse[1][1] / det_gesse;
				obr_gesse[0][1] = -1 * gesse[1][0] / det_gesse;
				obr_gesse[1][0] = -1 * gesse[0][1] / det_gesse;
				obr_gesse[1][1] = gesse[0][0] / det_gesse;
				cout << endl << "***�������� ������� �����***" << endl;
				for (int i = 0; i < 2; i++)
				{
					for (int j = 0; j < 2; j++)
						printf(" %7.3f", obr_gesse[i][j]);
					cout << endl;
				}
				//������������ �������� ������� ����� �� �������������������������
				if ((obr_gesse[0][0] > 0) && (det_gesse > 0))
				{
					cout << "�������� ������� ����� �����������������������" << endl;
					//���� ����������� ������
					ff[0] = dfx(x0, y0, 0, 0, 0);
					ff[1] = dfy(x0, y0, 0, 0, 0);
					for (int i = 0; i < 2; i++)
						for (int j = 0; j < 2; j++)
							dk[i] -= obr_gesse[i][j] * ff[j];
				}
				else
				{
					dk[0] = -dfx(x0, y0, 0, 0, 0);
					dk[1] = -dfy(x0, y0, 0, 0, 0);
				}
				cout << "------------------------------------" << endl;
				//���� �������� ���� s
				s0 = a;
				if (dfs(x0, y0, s0, dk[0], dk[1])*d3fs(x0, y0, s0, dk[0], dk[1]) > 0)  // ��� ������ ��������� ����� ��������� ������� �������� f(x0)*d2f(x0)>0
					s0 = a;
				else
					s0 = b;
				cout << "s0= " << s0 << endl;
				s = s0 - (f(x0, y0, s0, dk[0], dk[1]) / dfs(x0, y0, s0, dk[0], dk[1]));//������ �����������
				cout << "s= " << s << endl;
				while ((fabs(s0 - s) >= eps2)) // ���� �� ��������� ����������� ��������, ����� ���������� ���������
				{
					s0 = s;
					s = s0 - (dfs(x0, y0, s0, dk[0], dk[1]) / d2fs(x0, y0, s0, dk[0], dk[1])); //������� �������
				}

				//������� ����� ��������
				x = x0 + s * dk[0];
				y = y0 + s * dk[1];
				flag = 2;
			}
		} while (((fabs(f(x, y, 0, 0, 0) - f(x0, y0, 0, 0, 0)))>eps2) && ((sqrt((x - x0)*(x - x0) - (y - y0)*(y - y0)))>eps2));
		printf("\n->����� ��������: x = (%2.1f;%2.1f)\n", x, y);
		printf("\n->����� �������� ������� � ����� ��������: f(%2.1f;%2.1f) = %2.1f\n", x, y, f(x, y, 0, 0, 0));
		cout << endl;
		system("pause");
		return 0;
}


