#pragma once
#include <ilcplex/ilocplex.h>
#include "Problema.h"

class Modelo_Cplex :
	public Problema
{
	IloIntVarArray y_bi;
	IloArray<IloIntVarArray> y_tri;
	IloArray<IloBoolVarArray> x;

public:
	/*Usando o construtor do pai explicitamente*/
	Modelo_Cplex(const char * filename) : Problema(filename) {}

	~Modelo_Cplex();
};

