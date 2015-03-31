#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cassert>

class Matrix
{
	std::vector<double> data;
	int width;

public:
	double& operator()(int i, int j)
	{
		return data[(i - 1)*width + (j - 1)];
	}

	const double& operator()(int i, int j) const
	{
		return data[(i - 1)*width + (j - 1)];
	}

	double& operator[](int pos)
	{
		return data[(pos-1)];
	}

	const double& operator[](int pos) const
	{
		return data[(pos - 1)];
	}

	Matrix(int wiersze, int kolumny, double e = 0.0)
	{
		width = kolumny;
		data.resize(wiersze*kolumny, e);
	}

	int wiersze() const { return data.size() / width; }
	int kolumny() const { return width; }
};

void wypisz_macierz(const Matrix& m)
{
	for(int i = 1; i <= m.wiersze(); ++i)
	{
		for(int j = 1; j <= m.kolumny(); ++j)
		{
			std::cout << std::setw(10) << m(i, j) << " ";
		}
		std::cout << "\n";
	}
}

std::pair<Matrix, Matrix> dekompozycja_lu(const Matrix& a)
{
	assert(a.wiersze() == a.kolumny());
	int n = a.wiersze();
	std::pair<Matrix, Matrix> ret(Matrix(n, n), Matrix(n, n));
	Matrix& l = ret.first;
	Matrix& u = ret.second;

	// wype³nij jedynkami przek¹tn¹ macierzy L
	for(int i = 1; i <= n; ++i)
		l(i, i) = 1;

	for(int i = 1; i <= n; ++i)
	{
		for(int j = i; j <= n; ++j)
		{
			double suma = 0;
			for(int k = 1; k <= i - 1; ++k)
			{
				suma += l(i, k)*u(k, j);
			}
			u(i, j) = a(i, j) - suma;
		}
		for(int j = i + 1; j <= n; ++j)
		{
			double suma = 0;
			assert(u(i, i) != 0);
			for(int k = 1; k <= i - 1; ++k)
			{
				suma += l(j, k)*u(k, i);
			}
			l(j, i) = (a(j, i) - suma) / u(i, i);
		}
	}
	return ret;
}

Matrix rozwiaz_uklad(const Matrix& a, const Matrix& y)
{
	auto d = dekompozycja_lu(a);
	Matrix& l = d.first;
	Matrix& u = d.second;
	Matrix z(l.kolumny(), 1);
	Matrix x(l.kolumny(), 1);

	// l * z = y
	for(int i = 1; i <= l.wiersze(); ++i)
	{
		double suma = y[i];
		for(int j = 1; j <= i - 1; ++j)
			suma -= l(i, j)*z[j];
		z[i] = suma / l(i,i);
	}

	// u * x = z
	for(int i = l.wiersze(); i >= 1; --i)
	{
		double suma = z[i];
		for(int j = i + 1; j <= l.wiersze(); ++j)
			suma -= u(i, j)*x[j];
		x[i] = suma / u(i, i);
	}

	return x;
}

int main()
{
	const int op = 2;
	if(op == 0)
	{
		Matrix m(3, 2);
		m(1, 1) = 1;
		m(1, 2) = 2;
		m(2, 1) = 3;
		m(2, 2) = 4;
		m(3, 1) = 5;
		m(3, 2) = 6;
		wypisz_macierz(m);
	}
	else if(op == 1)
	{
		Matrix m(3, 3);
		m(1, 1) = 5;
		m(1, 2) = 3;
		m(1, 3) = 2;
		m(2, 1) = 1;
		m(2, 2) = 2;
		m(2, 3) = 0;
		m(3, 1) = 3;
		m(3, 2) = 0;
		m(3, 3) = 4;
		wypisz_macierz(m);
		auto d = dekompozycja_lu(m);
		wypisz_macierz(d.first);
		wypisz_macierz(d.second);
	}
	else if(op == 2)
	{
		Matrix m(3, 3);
		m(1, 1) = 5;
		m(1, 2) = 3;
		m(1, 3) = 2;
		m(2, 1) = 1;
		m(2, 2) = 2;
		m(2, 3) = 0;
		m(3, 1) = 3;
		m(3, 2) = 0;
		m(3, 3) = 4;

		Matrix y(1, 3);
		y[1] = 10;
		y[2] = 5;
		y[3] = -2;

		auto x = rozwiaz_uklad(m, y);
	}
}