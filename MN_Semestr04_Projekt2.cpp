#include <vector>
#include <utility>
#include <iostream>
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

	Matrix(int wiersze, int kolumny, double e = 0.0)
	{
		width = kolumny;
		data.resize(wiersze*kolumny, e);
	}

	int wiersze() const { return data.size() / width; }
	int kolumny() const { return width; }
};

int main()
{
	Matrix m(3, 2);
	m(1, 1) = 1;
	m(1, 2) = 2;
	m(2, 1) = 3;
	m(2, 2) = 4;
	m(3, 1) = 5;
	m(3, 2) = 6;
	for(int i = 1; i <= m.wiersze(); ++i)
	{
		for(int j = 1; j <= m.kolumny(); ++j)
		{
			std::cout << m(i, j) << " ";
		}
		std::cout << "\n";
	}
}