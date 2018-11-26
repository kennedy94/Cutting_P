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
	int	TamanhoDaPopulacao, NGeracoes, N_aleatorios;
	double	prob_mutacao, taxa_restart, taxa_elistimo;
	int crossover_media;
	default_random_engine generator;


protected:
	struct individuo {
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
		generator.seed(time(NULL));

		TamanhoDaPopulacao = 50;	//50, 100, 150
		NGeracoes = P * P; // P * P;		//100P, 300P, 500P, P^2
		prob_mutacao = 0.05;		//0.01, 0.05, 0.10
		taxa_elistimo = 0.50;		//0.10, 0.25, 0.50

		//Fixo depois dos primeiros taguchi
		taxa_restart = 0.2;
		N_aleatorios = 100 * P;
		crossover_media = 1; // 1;		
	}
	double fitness(individuo solu);

	individuo GerarSoluGRASP();

	bool viavel(individuo solu);

	void ImprimirVetorSolu(individuo solu);

	void selecao(vector<individuo>& Popu);

	individuo cruzar(individuo pai, individuo mae);

	individuo cruzar_diferenciado(individuo pai, individuo mae);

	void corrigir(individuo &solu);

	individuo cruzamento_unico(vector<individuo> Populacao);

	void mutar(individuo &solu);

	void Restart(vector<individuo> &Populacao);

	void ImprimirArquivo(individuo solu, double time);

	GA_Novo::individuo insert(individuo solu, int a, int b);

	individuo melhor_vizinho(individuo solu);

public:
	void Algoritmo_Genetico();

	void teste();

	GA_Novo(const char * filename) : Problema(filename) { };
};

