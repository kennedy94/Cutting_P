#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <list>
#include <limits>
#include <climits>
#include "Tipo_Viga.h"
#include "Padrao_Corte.h"
#include "Padrao_Pack.h"
#include "Padrao_Traspasse.h"

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
		P,			/*Número de padroes em empacotamento*/
		U,			/*Limite de leftovers*/
		O;			/*Número de padrões de traspasse*/
	double Maior_Forma, Menor_Forma;
	double Maior_Barra;
	double SPL_epsilon;
	int Maior_Valor_nos_Padroes;
	double Parameter_alpha_1, Parameter_alpha_2;
	vector<Tipo_Viga> TipoVigas;
	vector<Padrao_Corte> CutPatterns;
	vector<Padrao_Pack> PackPatterns;
	vector<Padrao_Traspasse> SplPatterns;
	vector<int> Gamma_Associado;

	vector<double> FORMAS, b, L, l, Menor_tamanho;
	vector<int> e;

	

public:
	Problema(const char * filename);
	double calculo_lowerbound();
	~Problema();
};

