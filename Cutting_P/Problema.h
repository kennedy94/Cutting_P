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
	int W,			/*N�mero de diferentes tamanhos de barras*/
		V,			/*N�mero de diferentes tamanhos de leftovers*/
		Gamma,		/*N�mero de diferentes tamanhos de formas*/
		M,			/*N�mero de formas dispon�veis*/
		C,			/*N�mero de tipos de vigas*/
		T,			/*N�mero de per�odos dispon�veis*/
		Maior_Qc,	/*Maior n�mero de tamanhos entre os tipos
					de viga*/
		H,			/*N�mero de padroes de corte*/
		P;			/*N�mero de padroes em empacotamento*/
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

