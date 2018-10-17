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
	double comparacao(individuo indie) {
		int contador = 0;

		for (int i = 0; i < this->ind.size(); i++)
			if (this->ind[i] == indie.ind[i] && this->n_vezes[i] == indie.n_vezes[i])
				contador++;

		return (double)contador / this->ind.size();
	}
	bool operator<(const individuo& rhs) const {
		return fitness < rhs.fitness;
	}
};

class Heuristica :
	public Problema
{
private:
	int	TamanhoDaPopulacao,
		NGeracoes;
	double	prob_mutacao,
			prob_cruzamento;
	double taxa_aumento_mut;
	vector<individuo> list_tabu;
			
protected:
	
	void definir_parametros() {
		TamanhoDaPopulacao = 30;
		prob_mutacao = 0.05;
		prob_cruzamento = 0.05;
		NGeracoes = 30;
		taxa_aumento_mut = 0.30;
	}


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

	void restart_populacao(vector<individuo>& Populacao);

	//testa se é uma solução viávels
	bool viavel(individuo solu);

	void ImprimirSolucaoEstiloCPLEX(individuo solu);

	//corrige uma solução inviável
	void corrigir(individuo &solu);

	int qtde_adicionavel(Padrao_Traspasse Padrao, vector<int> EstoqueUsado, list<int> usados);
	
	void ImprimirVetorSolu(individuo solu);

	void ImprimirArquivo(individuo solu, double time);

	individuo GerarSoluGRASP();

public:
	void funcaoteste();
	Heuristica(const char * filename) : Problema(filename) { };
};


