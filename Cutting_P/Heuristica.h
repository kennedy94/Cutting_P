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
			//taxa_aumento_mut;
			
protected:
	
	void definir_parametros() {
		TamanhoDaPopulacao = P;
		prob_mutacao = 0.05;
		//taxa_aumento_mut = 0.01;
		prob_cruzamento = 0.35;
		NGeracoes = 1000;
	}


	double fitness(individuo solu);

	//retorna v�rios indiv�duos obtidos por cruzando entre outros
	//cada um tem uma chance de muta��o
	list<individuo> cruzar(individuo pai, individuo mae);

	//cruza v�rios indiv�duos e retorna todos os gerados
	list<individuo> cruzamento(vector<individuo> Popu);

	/*junta os indiv�duos obtidos por cruzamento � popula��o atual e
	seleciona s� um n�mero predefinido deles
	*/
	void selecao(vector<individuo> &Popu);//

	//altera um elemento do vetor
	void mutar(individuo &solu);

	//testa se � uma solu��o vi�vel
	bool viavel(individuo solu);

	void ImprimirSolucaoEstiloCPLEX(individuo solu);

	//corrige uma solu��o invi�vel
	void corrigir(individuo &solu);

	int qtde_adicionavel(Padrao_Traspasse Padrao, vector<int> EstoqueUsado, list<int> usados);
	
	void ImprimirVetorSolu(individuo solu);

	void ImprimirArquivo(individuo solu, double time);

	individuo GerarSoluGRASP();

	individuo melhor_vizinho(individuo solu);

	individuo insert(individuo solu, int a, int b);

	void torneio(vector<individuo>& Popu);

	void Restart(vector<individuo>& Popu);

	

public:
	void funcaoteste();
	Heuristica(const char * filename) : Problema(filename) { };
};


