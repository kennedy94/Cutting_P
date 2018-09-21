#pragma once
#include "Problema.h"

class Heuristica :
	public Problema
{
private:
	list<individuo> Populacao;

protected:
	
	//aqui o bicho pega
	double fitness(individuo solu);

	//retorna vários indivíduos obtidos por cruzando entre outros
	//cada um tem uma chance de mutação
	list<individuo> cruzar(individuo pai, individuo mae);

	//cruza vários indivíduos e retorna todos os gerados
	list<individuo> cruzamento();

	/*junta os indivíduos obtidos por cruzamento à população atual e
	seleciona só um número predefinido deles
	*/
	void selecao();//

	//altera um elemento do vetor
	void mutar(individuo &solu);

	//testa se é uma solução viávels
	bool viavel(individuo solu);

	//corrige uma solução inviável
	void corrigir(individuo &solu);

public:
	Heuristica(const char * filename) : Problema(filename) { };
	~Heuristica();
};

struct individuo {
	vector<int> ind;
	vector<int> n_vezes;
};
