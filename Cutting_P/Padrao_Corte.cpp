#include "Padrao_Corte.h"



bool Padrao_Corte::Gera_LeftOver(int Gamma, int v)
{
	return tamanhos[Gamma + v] > 0; //lembra que t� zero based
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


