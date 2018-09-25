#pragma once
#include "Problema.h"
#include <numeric>
#include <random>

struct individuo {
	vector<int> ind;
	vector<int> n_vezes;
};

int kMenor(const vector<int> v, int k);

class Heuristica :
	public Problema
{
private:
	list<individuo> Populacao;

protected:
	
	//aqui o bicho pega
	double fitness(individuo solu);

	//retorna v�rios indiv�duos obtidos por cruzando entre outros
	//cada um tem uma chance de muta��o
	list<individuo> cruzar(individuo pai, individuo mae);

	//cruza v�rios indiv�duos e retorna todos os gerados
	list<individuo> cruzamento();

	/*junta os indiv�duos obtidos por cruzamento � popula��o atual e
	seleciona s� um n�mero predefinido deles
	*/
	void selecao();//

	//altera um elemento do vetor
	void mutar(individuo &solu);

	//testa se � uma solu��o vi�vels
	bool viavel(individuo solu);

	//corrige uma solu��o invi�vel
	void corrigir(individuo &solu);
	void ImprimirVetorSolu(individuo solu);
	individuo GerarSoluAleatoria();

public:
	void funcaoteste();
	Heuristica(const char * filename) : Problema(filename) { };
};


