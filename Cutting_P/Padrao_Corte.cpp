#include "Padrao_Corte.h"



bool Padrao_Corte::Gera_LeftOver(int W, int v)
{
	return tamanhos[W + v] > 0; //lembra que tá zero based
}

bool Padrao_Corte::Gera_LeftOvers(int Gamma, int V)
{
	for (int i = Gamma; i < Gamma + V; i++){
		if (tamanhos[i] == 0)
			continue;
		else
			return true;
	}
	return false;
}

Padrao_Corte::Padrao_Corte()
{
}


Padrao_Corte::~Padrao_Corte()
{
}


