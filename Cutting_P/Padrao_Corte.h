#pragma once
#include <vector>
#include <iostream>
using namespace std;

class Padrao_Corte
{
public:
	vector<int> tamanhos;
	double cap;
	int index_pat, index_barra;
	bool Gera_LeftOver(int W, int V);
	Padrao_Corte();
	~Padrao_Corte();
};

