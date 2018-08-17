#pragma once
#include <iostream>
#include <vector>
using namespace std;
class Tipo_Viga {
public:
	Tipo_Viga();
	~Tipo_Viga();
	const Tipo_Viga &operator=(const Tipo_Viga &obj);


	int tempo_cura, n_comprimentos, n_barras;
	vector<int> demandas;
	vector<double> comprimentos;
};
