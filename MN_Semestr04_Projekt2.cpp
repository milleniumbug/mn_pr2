#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>

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

Matrix& operator+=(Matrix& lhs, const Matrix& rhs)
{
	assert(lhs.wiersze() == rhs.wiersze());
	assert(lhs.kolumny() == rhs.kolumny());

	for(int i = 1; i <= lhs.wiersze(); ++i)
		for(int j = 1; j <= lhs.kolumny(); ++j)
			lhs(i, j) += rhs(i, j);
	return lhs;
}

Matrix operator+(Matrix lhs, const Matrix& rhs)
{
	lhs += rhs;
	return lhs;
}

Matrix& operator-=(Matrix& lhs, const Matrix& rhs)
{
	assert(lhs.wiersze() == rhs.wiersze());
	assert(lhs.kolumny() == rhs.kolumny());

	for(int i = 1; i <= lhs.wiersze(); ++i)
		for(int j = 1; j <= lhs.kolumny(); ++j)
			lhs(i, j) -= rhs(i, j);
	return lhs;
}

Matrix operator-(Matrix lhs, const Matrix& rhs)
{
	lhs -= rhs;
	return lhs;
}

int rozmiar(const Matrix& m)
{
	return m.wiersze() * m.kolumny();
}

Matrix macierz_skalarna(int n, double skalar)
{
	Matrix m(n, n);
	for(int i = 1; i <= n; ++i)
		m(i, i) = skalar;
	return m;
}

Matrix macierz_jednostkowa(int n)
{
	return macierz_skalarna(n, 1);
}

void swap_columns(Matrix& m, int a, int b)
{
	for(int i = 1; i <= m.wiersze(); ++i)
		std::swap(m(i, a), m(i, b));
}

void swap_rows(Matrix& m, int a, int b)
{
	for(int i = 1; i <= m.kolumny(); ++i)
		std::swap(m(a, i), m(b, i));
}

bool jest_diagonalna(const Matrix& m)
{
	for(int i = 1; i <= m.kolumny(); ++i)
		for(int j = 1; j <= m.wiersze(); ++j)
		{
			if(i == j)
				continue;
			if(m(i, j) != 0)
				return false;
		}
	return true;
}

// mno¿enie dwóch macierzy przy za³o¿eniu ¿e prawa jest diagonalna
Matrix matrix_multiplication_right_diagonal(Matrix x, const Matrix& b)
{
	assert(x.kolumny() == b.wiersze());
	assert(jest_diagonalna(b));
	//  amounts to multiplying the i-th column of A by ai.
	int n = std::min(b.wiersze(), b.kolumny());
	for(int i = 1; i <= n; ++i)
		for(int j = 1; j <= x.wiersze(); ++j)
			x(j, i) *= b(i, i);
	return x;
}

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
	assert([&]()
	{
		for(int i = 1; i <= a.wiersze(); ++i)
		{
			if(a(i, i) == 0)
				return false;
		}
		return true;
	}());
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

// "naprawia" macierze ¿eby mo¿na by³o je u¿yæ do obliczania uk³adów równañ razem z dekompozycj¹ LU
// i obliczaniem uk³adu równañ
Matrix fix_matrices(Matrix& a, Matrix& y)
{
	Matrix pozycje(rozmiar(y), 1);
	for(int i = 1; i <= y.wiersze(); ++i)
		pozycje[i] = i;

	for(int i = 1; i <= a.kolumny(); ++i)
	{
		int max_row = 1;
		for(int j = 1; j <= a.wiersze(); ++j)
		{
			if(std::abs(a(max_row, i)) <= std::abs(a(j, i)))
				max_row = j;
		}
		swap_rows(a, max_row, i);
		swap_rows(y, max_row, i);
		swap_rows(pozycje, max_row, i);
	}
	return pozycje;
}

Matrix rozwiaz_uklad(Matrix a, Matrix y)
{
	assert(a.wiersze() == a.kolumny());
	assert(a.wiersze() == y.wiersze());
	assert(y.kolumny() == 1);
	auto pozycje = fix_matrices(a, y);
	auto d = dekompozycja_lu(a);
	Matrix& l = d.first;
	Matrix& u = d.second;
	Matrix z(rozmiar(y), 1);
	

	// l * z = y
	for(int i = 1; i <= l.wiersze(); ++i)
	{
		double suma = y[i];
		for(int j = 1; j <= i - 1; ++j)
			suma -= l(i, j)*z[j];
		z[i] = suma / l(i,i);
	}

	Matrix x(rozmiar(y), 1);
	// u * x = z
	for(int i = l.wiersze(); i >= 1; --i)
	{
		double suma = z[i];
		for(int j = i + 1; j <= l.wiersze(); ++j)
			suma -= u(i, j)*x[j];
		x[i] = suma / u(i, i);
	}

	Matrix fixed_x(rozmiar(y), 1);
	for(int i = 1; i <= rozmiar(fixed_x); ++i)
	{
		fixed_x[pozycje[i]] = x[i];
	}
	return x;
}

int main()
{
	const int op = 6;
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

		Matrix y(3, 1);
		y[1] = 10;
		y[2] = 5;
		y[3] = -2;

		auto x = rozwiaz_uklad(m, y);
		wypisz_macierz(m);
		wypisz_macierz(y);
		wypisz_macierz(x);
	}
	else if(op == 3)
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

		swap_rows(m, 1, 3);
		swap_columns(m, 1, 3);
		wypisz_macierz(m);
	}
	else if(op == 4)
	{
		Matrix m(3, 3);
		m(1, 1) = 0;
		m(1, 2) = 0;
		m(1, 3) = 2;
		m(2, 1) = 7;
		m(2, 2) = 0;
		m(2, 3) = 0;
		m(3, 1) = 0;
		m(3, 2) = 4;
		m(3, 3) = 0;

		Matrix y(3, 1);
		y[1] = 10;
		y[2] = 5;
		y[3] = -2;

		auto poz = fix_matrices(m, y);
		wypisz_macierz(m);
		wypisz_macierz(y);
		wypisz_macierz(poz);
	}
	else if(op == 5)
	{
		Matrix m(3, 3);
		m(1, 1) = 0;
		m(1, 2) = 0;
		m(1, 3) = 2;
		m(2, 1) = 7;
		m(2, 2) = 0;
		m(2, 3) = 0;
		m(3, 1) = 0;
		m(3, 2) = 4;
		m(3, 3) = 0;

		Matrix n(3, 3);
		n(1, 1) = 1;
		n(1, 2) = 0;
		n(1, 3) = 0;
		n(2, 1) = 0;
		n(2, 2) = 4;
		n(2, 3) = 0;
		n(3, 1) = 0;
		n(3, 2) = 0;
		n(3, 3) = -1;

		std::cout << jest_diagonalna(m) << "\n";
		std::cout << jest_diagonalna(n) << "\n";
	}
	else if(op == 6)
	{
		
	}
}