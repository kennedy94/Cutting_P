#include "Modelo_Cplex.h"

void Modelo_Cplex::iniciar_variaveis(){
	Parameter_alpha_1 = 1;
	Parameter_alpha_2 = 1;

	//Definir os conjutos que geram ou não leftovers;
	for (IloInt h = 0; h < H; h++) {
		if (CutPatterns[h].Gera_LeftOvers(W, V))
			H_mais.push_back(h);
		else
			H_menos.push_back(h);
	}
	G.resize(Gamma);
	for (int gamma = 0; gamma < Gamma; gamma++) {
		for (int m = 0; m < M; m++)
			if (FORMAS[m] == L[gamma])
				G[gamma].push_back(m);
	}



	CPLEX_y_bi = IloArray<IloIntVarArray>(env, H_menos.size());
	for (auto h: H_menos)
		CPLEX_y_bi[h] = IloIntVarArray(env, W+V, 0, Maior_Valor_nos_Padroes);
	
	CPLEX_y_tri = IloArray<IloArray<IloIntVarArray>>(env, H_mais.size());
	for (auto h: H_mais) {
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


	for (auto h: H_menos)
		if (!CutPatterns[h].Gera_LeftOvers(W, V))
			for (IloInt w = 0; w < W; w++)
				soma1 += (b[w] - CutPatterns[h].cap) * CPLEX_y_bi[h][w];

	for (auto h: H_mais)
		for (IloInt w = 0; w < W; w++)
			for (IloInt v = 0; v < V; v++)
				if(CutPatterns[h].Gera_LeftOver(w,v))
					soma2 += (b[w] - CutPatterns[h].cap)
							* CPLEX_y_tri[h][w][v];


	for (auto h: H_menos)
		if (!CutPatterns[h].Gera_LeftOvers(W, V))
			for (IloInt w = W; w < W + V; w++)
				soma2 += (b[w] - CutPatterns[h].cap) * CPLEX_y_bi[h][w];

	model.add(IloMinimize(env, soma1 + Parameter_alpha_1*soma2 + 
		Parameter_alpha_2*soma3));
}


void Modelo_Cplex::restricoes_onlyone() {
	IloInt m, t, i;
	//para cada forma m e periodo de tempo t
	for (m = 0; m < M; m++) {
		for (t = 0; t < T; t++) {
			IloExpr expr(env);
			expr += CPLEX_x[0][m][t];
			for (i = 1; i < P; i++) {
				int tipo = PackPatterns[i].tipo;
				if ((PackPatterns[i].cap <= FORMAS[m]) 
					&& PackPatterns[i].maximal(FORMAS[m], TipoVigas[tipo].comprimentos)) {
					expr += CPLEX_x[i][m][t];
				}
			}
			model.add(expr <= 1).setName("Um Padrao");
			expr.end();
		}
	}
}


void Modelo_Cplex::restricoes_demanda() {
	IloInt c, k, m, t, i;
	//Para cada tipo de viga e tamanhos dentro do tipo de viga
	for (c = 0; c < C; c++) {
		for (k = 0; k < TipoVigas[c].n_comprimentos; k++) {
			//sum
			IloExpr expr(env);
			//problema: alguns problemas não estão gerando padrões com a última demanda
			for (m = 0; m < M; ++m) {
				for (i = 1; i < P; ++i) {
					int tipo = PackPatterns[i].tipo;
					for (t = 0; t < T - TipoVigas[PackPatterns[i].tipo].tempo_cura +1; ++t) {
						if (PackPatterns[i].cap <= FORMAS[m] && PackPatterns[i].tipo == c && PackPatterns[i].maximal(FORMAS[m], TipoVigas[tipo].comprimentos)) {
							expr += PackPatterns[i].tamanhos[k] * CPLEX_x[i][m][t];
						}
					}
				}
			}
			model.add(expr >= TipoVigas[c].demandas[k]).setName("Demanda");
			expr.end();
		}
	}
}

void Modelo_Cplex::restricoes_sequenciamento() {
	int i, j, m, a, t;
	for (m = 0; m < M; ++m)
		for (i = 1; i < P; ++i)
		{
			int tipo = PackPatterns[i].tipo;
			if (PackPatterns[i].cap <= FORMAS[m] && PackPatterns[i].maximal(FORMAS[m], TipoVigas[tipo].comprimentos))
				for (t = 0; t < T - TipoVigas[tipo].tempo_cura + 1; t++) {

					if (TipoVigas[tipo].tempo_cura != 1) {
						IloExpr expr(env);
						for (a = 1; a <= TipoVigas[tipo].tempo_cura - 1; a++)
							expr += CPLEX_x[0][m][t + a];
						model.add((TipoVigas[tipo].tempo_cura - 1) * CPLEX_x[i][m][t] <= expr).setName("problema");
						expr.end();
					}

				}
		}


	int R = 0;
	for (int c = 0; c < C; c++)
		if (TipoVigas[c].tempo_cura > R)
			R = TipoVigas[c].tempo_cura;

	for (t = 0; t < T; t++)
		for (m = 0; m < M; m++) {
			IloExpr expr(env);

			for (int beta = 2; beta <= R; beta++)
				for (j = beta; j <= R; j++)
					for (i = 0; i < P; i++) {
						int tipo = PackPatterns[i].tipo;
						if (PackPatterns[i].cap <= FORMAS[m] &&
							PackPatterns[i].maximal(FORMAS[m], TipoVigas[tipo].comprimentos) &&
							TipoVigas[PackPatterns[i].tipo].tempo_cura == j && t - beta + 1 >= 0)
						{
							expr += CPLEX_x[i][m][t - beta + 1];
						}

					}
			model.add(CPLEX_x[0][m][t] <= expr).setName("oi");
			expr.end();
		}
}

void Modelo_Cplex::restricoes_estoque() {
	//estoque de leftovers
	
	//leftovers cortados obedecem restricoes de estoque
	for (IloInt w = W; w < W + V; w++) {
		IloExpr expr(env);
		for (auto h: H_menos)
			expr += CPLEX_y_bi[h][w];

		model.add(expr <= e[w]);
		expr.end();
	}
	
	//barras novas cortadas obedecem o estoque
	for (IloInt w = 0; w < W; w++) {
		IloExpr expr1(env), expr2(env);
		for (auto h : H_menos)
			expr1 += CPLEX_y_bi[h][w];

		for (IloInt v = 0; v < V; w++)
			for (auto h : H_mais)
				expr1 += CPLEX_y_tri[h][w][v];

		model.add(expr1 + expr2 <= e[w]);
		expr1.end();
		expr2.end();
	}
}

//#medo
void Modelo_Cplex::restricoes_integracao() {

	for (IloInt gamma = 0; gamma < Gamma; gamma++)
	{
		IloExpr soma1(env), soma2(env), soma3(env), soma4(env);
		
		for (IloInt w = 0; w < W; w++)
			for (auto h : H_menos)
				if (CutPatterns[h].index_barra == w)
					soma1 += CutPatterns[h].tamanhos[gamma] * CPLEX_y_bi[h][w];

		for (IloInt w = 0; w < W; w++)
			for (IloInt v = 0; v < V; w++)
				for (auto h : H_mais)
					if (CutPatterns[h].index_barra == w 
						&& CutPatterns[h].Gera_LeftOver(w,v))
						soma2 +=
						CutPatterns[h].tamanhos[gamma] * CPLEX_y_tri[h][w][v];

		for (IloInt w = W; w < W+V; w++)
			for (auto h : H_menos)
				if(CutPatterns[h].index_barra == w)
					soma3 += CutPatterns[h].tamanhos[gamma] * CPLEX_y_bi[h][w];


		for (auto m: G[gamma])
			for (int t = 0; t < T; t++)
				for (int i = 1; i < P; i++)
					if (PackPatterns[i].cap <= FORMAS[m])
						soma4 += TipoVigas[PackPatterns[i].tipo].n_barras*CPLEX_x[i][m][t];

		model.add(soma1 + soma2 + soma3 == soma4);

		soma1.end();
		soma2.end();
		soma3.end();
		soma4.end();
	}

}
