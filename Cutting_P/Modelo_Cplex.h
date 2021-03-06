#pragma once
#include <ilcplex/ilocplex.h>
#include "Problema.h"
#include <chrono>
#include <ctime>

class Modelo_Cplex :
	public Problema
{
private:
	IloIntVarArray CPLEX_o;
	IloArray<IloIntVarArray> CPLEX_y_bi;
	IloArray<IloArray<IloIntVarArray>> CPLEX_y_tri;
	IloArray<IloArray<IloBoolVarArray>> CPLEX_x;
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
	IloBoolVarArray z;

	//IloNum Parameter_alpha_1, Parameter_alpha_2;

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

public:
	/*Usando o construtor do pai*/
	Modelo_Cplex(const char * filename) : Problema(filename) { };

	void MontarModelo();

	void ImprimirSolucaoArquivo();

	void ImprimirSolucao();

	void ImprimirGantt();

	void PlotarBarras();

	void CPLEX_objective_function_transformada();

	IloNum LB_makespan();

	double LB_custo();

	void resolver_inteira();

	~Modelo_Cplex();
};

