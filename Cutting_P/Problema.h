#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <list>
#include "Tipo_Viga.h"
#include "Padrao_Corte.h"
#include "Padrao_Pack.h"

using namespace std;

class Problema
{
protected:

	const char *nome_instancia;
	int W,			/*Número de diferentes tamanhos de barras*/
		V,			/*Número de diferentes tamanhos de leftovers*/
		Gamma,		/*Número de diferentes tamanhos de formas*/
		M,			/*Número de formas disponíveis*/
		C,			/*Número de tipos de vigas*/
		T,			/*Número de períodos disponíveis*/
		Maior_Qc,	/*Maior número de tamanhos entre os tipos
					de viga*/
		H,			/*Número de padroes de corte*/
		P;			/*Número de padroes em empacotamento*/
	double Maior_Forma, Menor_Forma;
	int Maior_Valor_nos_Padroes;
	vector<Tipo_Viga> TipoVigas;
	vector<Padrao_Corte> CutPatterns;
	vector<Padrao_Pack> PackPatterns;
	
	vector<double> FORMAS, b, L, l, Menor_tamanho;
	vector<int> e;

public:
	Problema(const char * filename);
	~Problema();
};

