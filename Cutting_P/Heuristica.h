#pragma once
#include "Problema.h"
#include <numeric>
#include <random>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <iomanip>


class Heuristica :
	public Problema
{
private:
	int	TamanhoDaPopulacao,
		NGeracoes;
	double	prob_mutacao,
		prob_cruzamento;
			//taxa_aumento_mut;

	default_random_engine generator;
			
protected:
	struct Utilizacao{
		int idx, util;
		bool operator<(const Utilizacao& rhs) const {
			return util < rhs.util;
		}
	};
	struct individuo {
		vector<int> ind;
		vector<int> n_vezes;
		double fitness;

		double distancia(individuo indie, int P, int H, int O) {
			int distancia = 0;
			vector<int>
				auxP(P - 1, 0),
				auxH(H, 0),
				auxO(O, 0);

			for (int i = 1; i < P; i++) {
				int valor_desse, valor_daquele;
				for (int j = 0; j < P - 1; j++)
					if (this->ind[j] == i) {
						valor_desse = this->n_vezes[j];
						break;
					}
				for (int j = 0; j < P - 1; j++)
					if (indie.ind[j] == i) {
						valor_daquele = indie.n_vezes[j];
						break;
					}
				distancia += abs(valor_desse - valor_daquele);
			}

			for (int i = 0; i < H; i++) {
				int valor_desse, valor_daquele;
				for (int j = P - 1; j < P - 1 + H; j++)
					if (this->ind[j] == i) {
						valor_desse = this->n_vezes[j];
						break;
					}
				for (int j = P - 1; j < P - 1 + H; j++)
					if (indie.ind[j] == i) {
						valor_daquele = indie.n_vezes[j];
						break;
					}
				distancia += abs(valor_desse - valor_daquele);
			}



			for (int i = 0; i < O; i++) {
				int valor_desse, valor_daquele;
				for (int j = P - 1 + H; j < P - 1 + H + O; j++)
					if (this->ind[j] == i) {
						valor_desse = this->n_vezes[j];
						break;
					}
				for (int j = P - 1 + H; j < P - 1 + H + O; j++)
					if (indie.ind[j] == i) {
						valor_daquele = indie.n_vezes[j];
						break;
					}

				distancia += abs(valor_desse - valor_daquele);
			}

			return distancia;
		}

		bool operator<(const individuo& rhs) const {
			return fitness < rhs.fitness;
		}
	};
	void definir_parametros() {
		TamanhoDaPopulacao = 50;
		prob_mutacao = 0.05;
		prob_cruzamento = 0.25;
		NGeracoes = 5000;
		generator.seed(time(NULL));
		srand(time(NULL));
	}

	double distancia_total_popu(vector<individuo> Popu);
	double fitness(individuo solu);

	//retorna vários indivíduos obtidos por cruzando entre outros
	//cada um tem uma chance de mutação
	list<individuo> cruzar(individuo pai, individuo mae);

	//cruza vários indivíduos e retorna todos os gerados
	list<individuo> cruzamento(vector<individuo> Popu);

	list<individuo> cruzamento_diferenciado(vector<individuo> Popu);

	/*junta os indivíduos obtidos por cruzamento à população atual e
	seleciona só um número predefinido deles
	*/
	void selecao(vector<individuo> &Popu);//

	//altera um elemento do vetor
	void mutar(individuo &solu);

	//testa se é uma solução viável
	bool viavel(individuo solu);

	void ImprimirSolucaoEstiloCPLEX(individuo solu);

	//corrige uma solução inviável
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


