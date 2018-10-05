#pragma once
#include "Problema.h"
#include <numeric>
#include <random>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <iomanip>
struct individuo {
	vector<int> ind;
	vector<int> n_vezes;
	double fitness;
};

inline int kMenor(const vector<int> v, int k);

inline int kMaior(const vector<int> v, int k);

class Heuristica :
	public Problema
{
private:
	list<individuo> Populacao;
	int TamanhoDaPopulacao;
protected:
	
	double fitness(individuo solu);

	//retorna vários indivíduos obtidos por cruzando entre outros
	//cada um tem uma chance de mutação
	list<individuo> cruzar(individuo pai, individuo mae);

	//cruza vários indivíduos e retorna todos os gerados
	list<individuo> cruzamento(vector<individuo> Popu);

	/*junta os indivíduos obtidos por cruzamento à população atual e
	seleciona só um número predefinido deles
	*/
	void selecao(vector<individuo> &Popu);//

	//altera um elemento do vetor
	void mutar(individuo &solu);

	//testa se é uma solução viávels
	bool viavel(individuo solu);

	//corrige uma solução inviável
	void corrigir(individuo &solu);
	int qtde_adicionavel(Padrao_Traspasse Padrao, vector<int> EstoqueUsado, list<int> usados);
	void ImprimirVetorSolu(individuo solu);
	individuo GerarSoluGRASP();

public:
	void funcaoteste();
	Heuristica(const char * filename) : Problema(filename) { };
};


