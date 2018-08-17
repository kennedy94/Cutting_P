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

	IloEnv env;
	IloModel model;

	IloNum Parameter_alpha_1, Parameter_alpha_2;

	void CPLEX_objective_function();

public:
	/*Usando o construtor do pai*/
	Modelo_Cplex(const char * filename) : Problema(filename) {}

	void iniciar_variaveis();

	~Modelo_Cplex();
};

