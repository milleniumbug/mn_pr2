#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <string>

class Matrix
{
	typedef std::vector<double> underlying_container;
	underlying_container data;
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
		return data[(pos - 1)];
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

	for(int i = 1; i <= lhs.kolumny(); ++i)
		for(int j = 1; j <= lhs.wiersze(); ++j)
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

	for(int i = 1; i <= lhs.kolumny(); ++i)
		for(int j = 1; j <= lhs.wiersze(); ++j)
			lhs(i, j) -= rhs(i, j);
	return lhs;
}

Matrix operator-(Matrix lhs, const Matrix& rhs)
{
	lhs -= rhs;
	return lhs;
}

Matrix& operator*=(Matrix& lhs, double x)
{
	for(int i = 1; i <= lhs.kolumny(); ++i)
		for(int j = 1; j <= lhs.wiersze(); ++j)
			lhs(i, j) *= x;
	return lhs;
}

Matrix operator*(Matrix m, double x)
{
	m *= x;
	return m;
}

Matrix operator*(double x, Matrix m)
{
	m *= x;
	return m;
}

Matrix mnozenie_przez_skalar(Matrix a, double x)
{
	for(int i = 1; i <= a.kolumny(); ++i)
		for(int j = 1; j <= a.wiersze(); ++j)
			a(i, j) *= x;
	return a;
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

void zamien_kolumny(Matrix& m, int a, int b)
{
	for(int i = 1; i <= m.wiersze(); ++i)
		std::swap(m(i, a), m(i, b));
}

void zamien_rzedy(Matrix& m, int a, int b)
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

bool jest_kwadratowa(const Matrix& a)
{
	return a.wiersze() == a.kolumny();
}

// mno�enie dw�ch macierzy przy za�o�eniu �e prawa jest diagonalna
void mnozenie_prawa_diagonalna_inplace(Matrix& x, const Matrix& b)
{
	assert(jest_kwadratowa(x));
	assert(jest_diagonalna(b));
	//  amounts to multiplying the i-th column of A by ai.
	int n = std::min(b.wiersze(), b.kolumny());
	for(int i = 1; i <= n; ++i)
		for(int j = 1; j <= n; ++j)
			x(j, i) *= b(i, i);
}

Matrix mnozenie_prawa_diagonalna(Matrix x, const Matrix& b)
{
	mnozenie_prawa_diagonalna_inplace(x, b);
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
	assert(jest_kwadratowa(a));
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

	// wype�nij jedynkami przek�tn� macierzy L
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

// "naprawia" macierze �eby mo�na by�o je u�y� do obliczania uk�ad�w r�wna� razem z dekompozycj� LU
// i obliczaniem uk�adu r�wna�
Matrix napraw_macierze(Matrix& a, Matrix& y)
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
		zamien_rzedy(a, max_row, i);
		zamien_rzedy(y, max_row, i);
		zamien_rzedy(pozycje, max_row, i);
	}
	return pozycje;
}

Matrix rozwiaz_uklad(Matrix a, Matrix y)
{
	assert(jest_kwadratowa(a));
	assert(a.wiersze() == y.wiersze());
	assert(y.kolumny() == 1);
	auto pozycje = napraw_macierze(a, y);
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
		z[i] = suma / l(i, i);
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

Matrix page_rank(Matrix B, const double d)
{
	assert(jest_kwadratowa(B));
	const int n = B.wiersze();
	Matrix A(n, n);
	Matrix W(n, 1);
	for(int i = 1; i <= n; ++i)
	{
		int ilosc_linkow = 0;
		for(int j = 1; j <= n; ++j)
		{
			ilosc_linkow += B(j, i) ? 1 : 0;
		}
		if(ilosc_linkow == 0)
		{
			for(int j = 1; j <= n; ++j)
			{
				if(j == i)
					continue;
				else
					B(j, i) = 1;
			}
			ilosc_linkow = n - 1;
		}
		A(i, i) = 1.0 / ilosc_linkow;
	}

	for(int i = 1; i <= n; ++i)
	{
		W[i] = (1 - d) / n;
	}

	Matrix L = macierz_jednostkowa(n);
	mnozenie_prawa_diagonalna_inplace(B, A);
	B *= d;
	L -= B;
	return rozwiaz_uklad(L, W);
}

void testy(int op)
{
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

		zamien_rzedy(m, 1, 3);
		zamien_kolumny(m, 1, 3);
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

		auto poz = napraw_macierze(m, y);
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
		const int q = 1;
		const int w = 2;
		const int e = 3;
		const int r = 4;

		const double d = 0.85;
		Matrix B(4, 4);
		// Za��my przyk�adowo, �e ca�a sie� sk�ada si� z czterech stron o nazwach q, w, e i r.Niech strona
		// w ma link do strony q i strony e, strona e ma link do strony q, i strona r posiada linki do stron
		// q, w i e
		B(q, w) = 1;
		B(e, w) = 1;
		B(q, e) = 1;
		B(q, r) = 1;
		B(w, r) = 1;
		B(e, r) = 1;

		auto R = page_rank(B, d);
		wypisz_macierz(R);
	}
}

struct GraphReader
{
	std::istream& f;
	bool operator()(std::pair<int, int>& edge)
	{
		while(f.peek() == '#')
			f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		f >> edge.first;
		if(f.eof())
			return false;
		assert(!f.fail() && !f.bad());
		f >> edge.second;
		if(f.eof())
			return false;
		assert(!f.fail() && !f.bad());
		return true;
	}

	explicit GraphReader(std::istream& i) : f(i) {}
};

int najwiekszy_wierzcholek_w_pliku_z_grafem(std::istream& f)
{
	int max_vertex = -1;
	GraphReader r(f);

	std::pair<int, int> edge;
	while(r(edge))
	{
		max_vertex = std::max(max_vertex, edge.first);
		max_vertex = std::max(max_vertex, edge.second);
	}
	return max_vertex;
}

template<typename T>
void wczytaj_i_oblicz_page_rank(T f)
{
	std::string path;
	std::cin.clear();
	std::cin.ignore(256, '\n');
	std::cout << "Podaj sciezke do pliku\n";
	std::getline(std::cin, path);
	int max_vertex;
	{
		std::ifstream plik(path);
		max_vertex = najwiekszy_wierzcholek_w_pliku_z_grafem(plik);
	}
	Matrix m(max_vertex, max_vertex);
	{
		std::ifstream plik(path);
		GraphReader r(plik);

		std::pair<int, int> edge;
		while(r(edge))
		{
			f(edge, m);
		}
	}
	wypisz_macierz(page_rank(m, 0.85));
}

int main()
{
	std::cout
		<< "1. Obliczanie ukladow rownan\n"
		<< "2. Obliczanie PageRank dla grafow skierowanych\n"
		<< "3. Obliczanie PageRank dla grafow nieskierowanych\n";
	int select;
	std::cin >> select;
	if(select == 1)
	{
		int ilosc_rownan;
		std::cout << "Podaj ilosc rownan\n";
		std::cin >> ilosc_rownan;
		std::cout << "Uklad rownan A*x = y, wyznaczanie x.\n";
		std::cout << "Podaj wspolczynniki dla macierzy:\n";
		std::cout
			<< "( na przyklad (dla ilosci rownan 3):\n"
			<< " 1 2 3\n"
			<< " 4 5 6\n"
			<< " 7 8 9 )\n\n";
		Matrix m(ilosc_rownan, ilosc_rownan);
		Matrix y(ilosc_rownan, 1);
		for(int i = 1; i <= ilosc_rownan; ++i)
			for(int j = 1; j <= ilosc_rownan; ++j)
				std::cin >> m(i, j);
		std::cout << "Podaj y:\n";
		for(int i = 1; i <= ilosc_rownan; ++i)
			std::cin >> y[i];

		std::cout << "Obliczone x:\n";
		wypisz_macierz(rozwiaz_uklad(m, y));
	}
	else if(select == 2)
	{
		wczytaj_i_oblicz_page_rank([](std::pair<int, int> edge, Matrix& m)
		{
			m(edge.second, edge.first) = 1;
		});
	}
	else if(select == 3)
	{
		wczytaj_i_oblicz_page_rank([](std::pair<int, int> edge, Matrix& m)
		{
			m(edge.second, edge.first) = 1;
			m(edge.first, edge.second) = 1;
		});
	}
	std::cin.get();
}