#pragma once
#include <ilcplex/ilocplex.h>
#include "Problema.h"


class Modelo_Cplex :
	public Problema
{
private:
	IloArray<IloIntVarArray> CPLEX_y_bi;
	IloArray<IloArray<IloIntVarArray>> CPLEX_y_tri;
	IloArray<IloArray<IloBoolVarArray>> CPLEX_x;
	list<int> 
		H_mais,	//Guarda os indices dos padroes de
						//corte que geram leftover
		H_menos;//Guarda os indices dos padroes de
						//corte que nao geram leftover
	vector<list<int>> G; /*Vetor de conjuntos que guardam as formas
						 de tamanho gamma*/
		

	IloEnv env;
	IloModel model;

	IloNum Parameter_alpha_1, Parameter_alpha_2;

	void CPLEX_objective_function();

	void restricoes_onlyone();

	void restricoes_demanda();


	void restricoes_sequenciamento();

	void restricoes_estoque();

	void restricoes_integracao();

public:
	/*Usando o construtor do pai*/
	Modelo_Cplex(const char * filename) : Problema(filename) {}

	void iniciar_variaveis();

	~Modelo_Cplex();
};

