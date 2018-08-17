#include "Padrao_Pack.h"

Padrao_Pack::Padrao_Pack() {
	id = -1;
	tipo = -1;
	n_cobertos = 0;
	n_comprimentos = 0;
	cap = 0;
}

void Padrao_Pack::contar() {
	n_cobertos = 0;
	for (int i = 0; i < n_comprimentos; i++)
		if (tamanhos[i] > 0) { n_cobertos++; }
}
//Calcula o numero de tamanhos cobertos pelo padrao
void Padrao_Pack::gerar_cobertos(int k) {
	n_cobertos = 0;
	for (int i = 0; i < k; i++)
		if (tamanhos[i] > 0) n_cobertos++;

}
//Recebe o indice do tamanho e retorna se o Padrao_Pack o contem
bool Padrao_Pack::contem(int tam) {
	return (tamanhos[tam] > 0);
}


void Padrao_Pack::alocar_PADRAO(int k, int tipo) {
	this->n_comprimentos = k;
	this->tipo = tipo;
	tamanhos.resize(k);
	for (int i = 0; i < k; i++)
		tamanhos[i] = 0;

}

bool Padrao_Pack::comparar_demandas(const Tipo_Viga & c1)
{
	if (c1.n_comprimentos != n_comprimentos)	return false;

	for (int i = 0; i < c1.n_comprimentos; i++) {
		if (c1.demandas[i] > tamanhos[i])
			return false;
		else
			continue;
	}
	return true;
}

bool Padrao_Pack::comparar_demandas(const Tipo_Viga & c1, int IND_TAMANHO)
{
	if (c1.n_comprimentos != n_comprimentos)	return false;

	if (c1.demandas[IND_TAMANHO] > tamanhos[IND_TAMANHO])
		return false;

	return true;
}


bool operator== (const Padrao_Pack &c1, const Padrao_Pack &c2) {
	if (c1.tipo != c2.tipo) return false;

	for (int i = 0; i < c1.n_comprimentos; i++)
		if (c1.tamanhos[i] != c2.tamanhos[i]) return false;

	return true;
}
bool operator< (const Padrao_Pack &c1, const Padrao_Pack &c2) {
	return (c1.n_cobertos < c2.n_cobertos);
}

bool operator> (const Padrao_Pack &c1, const Padrao_Pack &c2) {
	return (c1.n_cobertos > c2.n_cobertos);
}

bool operador_padrao(const Padrao_Pack &c1, const Padrao_Pack &c2) {
	return (c1.n_cobertos > c2.n_cobertos);
}

bool operador_padrao_naocobertos(const Padrao_Pack &c1, const Padrao_Pack &c2) {
	return (c1.n_cobre_naocobertos > c2.n_cobre_naocobertos);
}

