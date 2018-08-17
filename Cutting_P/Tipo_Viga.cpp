#include "Tipo_Viga.h"



Tipo_Viga::Tipo_Viga()
{
}


Tipo_Viga::~Tipo_Viga()
{
}

const Tipo_Viga & Tipo_Viga::operator=(const Tipo_Viga &obj) {
	if (this == &obj) return *this;

	demandas.clear();
	comprimentos.clear();

	tempo_cura = obj.tempo_cura;
	n_comprimentos = obj.n_comprimentos;

	demandas.resize(n_comprimentos);
	comprimentos.resize(n_comprimentos);

	for (int i = 0; i < n_comprimentos; i++){
		demandas[i] = obj.demandas[i];
		comprimentos[i] = obj.comprimentos[i];
	}

	return *this;
}