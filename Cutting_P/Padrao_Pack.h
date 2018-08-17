#pragma once
#include <iostream>
#include "Tipo_Viga.h"
using namespace std;



class Padrao_Pack {
protected:
	void gerar_cobertos(int k);


public:

	//Construtores
	Padrao_Pack();
	//Estruturasaa
	int id, 
		tipo,
		n_cobertos,
		n_comprimentos;
	double cap;
	vector<int> tamanhos;
	
	//Métodos

	bool comparar_demandas(const Tipo_Viga & c1);
	bool comparar_demandas(const Tipo_Viga & c1, int IND_TAMANHO);
	void contar();
	bool contem(int tam);
	int n_cobre_naocobertos;
	void alocar_PADRAO(int k, int tipo);


};

bool operator== (const Padrao_Pack &c1, const Padrao_Pack &c2);

bool operator< (const Padrao_Pack &c1, const Padrao_Pack &c2);

bool operator> (const Padrao_Pack &c1, const Padrao_Pack &c2);

bool operador_padrao(const Padrao_Pack &c1, const Padrao_Pack &c2);

bool operador_padrao_naocobertos(const Padrao_Pack & c1, const Padrao_Pack & c2);