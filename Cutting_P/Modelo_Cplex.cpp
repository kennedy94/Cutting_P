#include "Modelo_Cplex.h"

void Modelo_Cplex::iniciar_variaveis(){
	Parameter_alpha_1 = 1;
	Parameter_alpha_2 = 1;
	CPLEX_y_bi = IloArray<IloIntVarArray>(env, H);
	for (IloInt h = 0; h < H; h++) 
		CPLEX_y_bi[h] = IloIntVarArray(env, W, 0, Maior_Valor_nos_Padroes);
	
	CPLEX_y_tri = IloArray<IloArray<IloIntVarArray>>(env, H);
	for (IloInt h = 0; h < H; h++) {
		CPLEX_y_tri[h] = IloArray<IloIntVarArray>(env, W);
		for (IloInt w = 0; w < W; w++)
			CPLEX_y_tri[h][w] = IloIntVarArray(env, V, 0, Maior_Valor_nos_Padroes);
	}
	
	CPLEX_x = IloArray<IloArray<IloBoolVarArray>>(env, P);
	for (IloInt i = 0; i < P; i++) {
		CPLEX_x[i] = IloArray<IloBoolVarArray>(env, M);
		for (IloInt m = 0; m < M; m++)
			CPLEX_x[i][m] = IloBoolVarArray(env, T);
	}
}

Modelo_Cplex::~Modelo_Cplex()
{
}

void Modelo_Cplex::CPLEX_objective_function(){
	IloExpr soma1, soma2, soma3;

	for (IloInt h = 0; h < H; h++)
		if (!Padroes_de_Corte[h].Gera_LeftOvers(W, V))
			for (IloInt w = 0; w < W; w++)
				soma1 += (b[w] - Padroes_de_Corte[h].cap) * CPLEX_y_bi[h][w];

	for (IloInt h = 0; h < H; h++)
		for (IloInt w = 0; w < W; w++)
			for (IloInt v = 0; v < V; v++)
				if(Padroes_de_Corte[h].Gera_LeftOver(w,v))
					soma2 += (b[w] - Padroes_de_Corte[h].cap)
							* CPLEX_y_tri[h][w][v];


	for (IloInt h = 0; h < H; h++)
		if (!Padroes_de_Corte[h].Gera_LeftOvers(W, V))
			for (IloInt w = W; w < W + V; w++)
				soma2 += (b[w] - Padroes_de_Corte[h].cap) * CPLEX_y_bi[h][w];



	model.add(IloMinimize(env, soma1 + Parameter_alpha_1*soma2 + 
		Parameter_alpha_2*soma3));
}

