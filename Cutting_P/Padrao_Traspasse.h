#pragma once
#include <iostream>
#include <vector>

using namespace std;

class Padrao_Traspasse
{
public:
	int id,
		barra_gerada;
	double folga;
	vector<int> tamanhos;
	bool operator<(const Padrao_Traspasse& rhs) const {
		return folga < rhs.folga;
	}


	Padrao_Traspasse();
	~Padrao_Traspasse();
};

