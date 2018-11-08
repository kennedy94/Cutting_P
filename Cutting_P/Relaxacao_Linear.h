#pragma once
#include <ilcplex/ilocplex.h>
#include "Problema.h"
#include <chrono>
#include <ctime>

class Relaxacao_Linear :
	public Problema
{
private:
	inline void CPLEX_objective_function();

	inline void restricoes_onlyone();

	inline void restricoes_demanda();

	inline void restricoes_z();

	inline void restricoes_sequenciamento();

	inline void restricoes_estoque();

	inline void restricoes_continuidade();

	inline void restricoes_integracao();

	void restricao_limite();

	void iniciar_variaveis();

	IloNumVarArray CPLEX_o;
	IloArray<IloNumVarArray> CPLEX_y_bi;
	IloArray<IloArray<IloNumVarArray>> CPLEX_y_tri;
	IloArray<IloArray<IloNumVarArray>> CPLEX_x;
	IloNumVar Makespan, Custo;


	list<int>
		H_mais,	//Guarda os indices dos padroes de
				//corte que geram leftover
		H_menos;//Guarda os indices dos padroes de
				//corte que nao geram leftover
	vector<list<int>> G; /*Vetor de conjuntos que guardam as formas
						 de tamanho gamma*/

	IloEnv env;
	double tempo_solucao;
	IloModel model;
	IloCplex cplex;
	IloNumVarArray z;
public:
	Relaxacao_Linear(const char * filename) : Problema(filename) { };

	void MontarModelo();

	void ImprimirSolucaoArquivo();

	double calculo_lowerbound();

	void resolver_linear();

};





