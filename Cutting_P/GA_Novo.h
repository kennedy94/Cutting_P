#pragma once
#include "Problema.h"
#include <numeric>
#include <random>
#include <chrono>
#include <time.h>

class GA_Novo :
	public Problema
{
private:
	int	TamanhoDaPopulacao,NGeracoes;
	double	prob_mutacao;
	default_random_engine generator;


protected:
	struct individuo{
		vector<int> ind;
		vector<int> n_vezes;
		double fitness;

		bool operator<(const individuo& rhs) const {
			return fitness < rhs.fitness;
		}


		bool igual(individuo indie) {
			if (indie.ind.size() != this->ind.size())
				return false;
			
			for (int i = 0; i < indie.ind.size(); i++)
				if (!(indie.ind[i] == this->ind[i] && indie.n_vezes[i] == this->n_vezes[i]))
					return false;

			return true;
		}


	};
	struct Utilizacao {
		int idx, util;
		bool operator<(const Utilizacao& rhs) const {
			return util < rhs.util;
		}
	};
	
	void definir_parametros() {
		srand(time(NULL));
		TamanhoDaPopulacao = 50;
		prob_mutacao = 0.1;
		NGeracoes = 100*P;
		generator.seed(time(NULL));
	}
	double fitness(individuo solu);

	individuo GerarSoluGRASP();

	bool viavel(individuo solu);

	void ImprimirVetorSolu(individuo solu);

	void selecao(vector<individuo>& Popu);

	individuo cruzar(individuo pai, individuo mae);

	void corrigir(individuo &solu);

	individuo cruzamento_unico(vector<individuo> Populacao);

	void mutar(individuo &solu);

	void Restart(vector<individuo> &Populacao);

	void ImprimirArquivo(individuo solu, double time);

	individuo melhor_vizinho(individuo solu);

public:
	void funcao_teste();
	
	GA_Novo(const char * filename) : Problema(filename) { };
};

