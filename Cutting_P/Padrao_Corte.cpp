#include "Padrao_Corte.h"



bool Padrao_Corte::Gera_LeftOver(int w, int v)
{
	return tamanhos[w + v - 1] > 0; //lembra que tá zero based
}

Padrao_Corte::Padrao_Corte()
{
}


Padrao_Corte::~Padrao_Corte()
{
}


