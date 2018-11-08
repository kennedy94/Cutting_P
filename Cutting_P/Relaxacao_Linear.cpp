#include "Relaxacao_Linear.h"


void Relaxacao_Linear::iniciar_variaveis() {

	char strnum[30];

	z = IloNumVarArray(env, T, 0, 1);
	for (int t = 0; t < T; t++) {
		sprintf(strnum, "z(%d)", t);
		z[t].setName(strnum);
	}


	CPLEX_o = IloNumVarArray(env, O, 0, 10000000);

	for (int mu = 0; mu < O; mu++) {
		sprintf(strnum, "o(%d)", mu);
		CPLEX_o[mu].setName(strnum);
	}

	//Definir os conjutos que geram ou não leftovers;
	for (IloInt h = 0; h < H; h++) {
		if (CutPatterns[h].cap > b[CutPatterns[h].index_barra]) {
			cout << CutPatterns[h].index_barra << " :";
			for (int i = 0; i < Gamma + V; i++)
				cout << CutPatterns[h].tamanhos[i] << " ";
			getchar();
			cout << endl;
		}

		if (CutPatterns[h].Gera_LeftOvers(Gamma, V)) {
			H_mais.push_back(h);
		}
		else {
			H_menos.push_back(h);

		}
	}

	G.resize(Gamma);
	for (int gamma = 0; gamma < Gamma; gamma++) {
		for (int m = 0; m < M; m++)
			if (FORMAS[m] == L[gamma])
				G[gamma].push_back(m);
	}

	CPLEX_y_bi = IloArray<IloNumVarArray>(env, H);
	for (int h = 0; h < H; h++) {
		CPLEX_y_bi[h] = IloNumVarArray(env, W + V, 0, 10000000);
		for (IloInt w = 0; w < W + V; w++)
		{
			sprintf(strnum, "y(%d,%d)", h, w);
			CPLEX_y_bi[h][w].setName(strnum);
		}
	}



	CPLEX_y_tri = IloArray<IloArray<IloNumVarArray>>(env, H);
	for (int h = 0; h < H; h++) {
		CPLEX_y_tri[h] = IloArray<IloNumVarArray>(env, W);
		for (IloInt w = 0; w < W; w++) {
			CPLEX_y_tri[h][w] = IloNumVarArray(env, V, 0, 10000000);
			for (IloInt v = 0; v < V; v++)
			{
				sprintf(strnum, "y(%d,%d,%d)", h, w, v);
				CPLEX_y_tri[h][w][v].setName(strnum);
			}
		}
	}

	CPLEX_x = IloArray<IloArray<IloNumVarArray>>(env, P);
	for (IloInt i = 0; i < P; i++) {
		CPLEX_x[i] = IloArray<IloNumVarArray>(env, M);
		for (IloInt m = 0; m < M; m++) {
			CPLEX_x[i][m] = IloNumVarArray(env, T, 0, 1);
			for (IloInt t = 0; t < T; t++)
			{
				sprintf(strnum, "x(%d,%d,%d)", i, m, t);
				CPLEX_x[i][m][t].setName(strnum);
			}
		}
	}
}

void Relaxacao_Linear::resolver_linear() {
	cplex = IloCplex(model);

	try
	{
		//cplex.setParam(IloCplex::PreInd, 0);
		//cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 8000);
		
		//cplex.setParam(IloCplex::NumericalEmphasis, 1);
		cplex.setParam(IloCplex::TiLim, 3600);
		auto start = chrono::system_clock::now();
		cplex.solve();
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		tempo_solucao = elapsed_seconds.count();

	}
	catch (IloException& e) {
		cerr << endl << e.getMessage() << endl;
		system("pause");
		throw(0);
	}
	catch (exception& e)
	{
		cerr << endl << e.what() << endl;
		system("pause");
		throw(1);
	}
	catch (...) {
		cerr << endl << "Socorro" << endl;
		system("pause");
	}

}

void Relaxacao_Linear::CPLEX_objective_function() {
	IloExpr soma1(env), soma2(env), soma3(env), expr(env);


	IloInt t;

	IloExpr costSum(env);
	for (t = 0; t < T; t++)
		costSum += z[t];


	for (auto h : H_menos) {
		if (CutPatterns[h].index_barra < W) {
			soma1 += (b[CutPatterns[h].index_barra] - CutPatterns[h].cap) * CPLEX_y_bi[h][CutPatterns[h].index_barra];
		}
	}
	for (auto mu : SplPatterns) {
		expr += mu.folga * CPLEX_o[mu.id];
	}

	for (auto h : H_mais) {
		if (CutPatterns[h].index_barra < W) {
			for (IloInt v = 0; v < V; v++) {
				if (CutPatterns[h].Gera_LeftOver(Gamma, v)) {
					soma2 += (b[CutPatterns[h].index_barra] - CutPatterns[h].cap) * CPLEX_y_tri[h][CutPatterns[h].index_barra][v];
				}
			}
		}
	}

	for (auto h : H_menos) {
		if (CutPatterns[h].index_barra >= W) {
			soma3 += (b[CutPatterns[h].index_barra] - CutPatterns[h].cap) * CPLEX_y_bi[h][CutPatterns[h].index_barra];
		}
	}
	//costSum +
	Makespan = IloNumVar(env, 0, T);
	Custo = IloNumVar(env, 0, 99999999);

	model.add(costSum == Makespan);
	model.add(Custo == soma1 + Parameter_alpha_1*soma2 + Parameter_alpha_2*(soma3 + expr));
	model.add(IloMinimize(env, Makespan + Custo));

	soma1.end();
	soma2.end();
	soma3.end();
	expr.end();
	costSum.end();
}

//Inteira
void Relaxacao_Linear::restricoes_onlyone() {
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



//SOS
//void Relaxacao_Linear::restricoes_onlyone() {
//	IloInt m, t, i;
//	//para cada forma m e periodo de tempo t
//	for (m = 0; m < M; m++) {
//		for (t = 0; t < T; t++) {
//			IloExpr expr(env);
//			IloNumVarArray teste(env);
//			teste.add(CPLEX_x[0][m][t]);
//			for (i = 1; i < P; i++) {
//				int tipo = PackPatterns[i].tipo;
//				if ((PackPatterns[i].cap <= FORMAS[m])
//					&& PackPatterns[i].maximal(FORMAS[m], TipoVigas[tipo].comprimentos)) {
//					teste.add(CPLEX_x[i][m][t]);
//				}
//			}
//			model.add(IloSOS1(env, teste));
//			//model.add(expr <= 1).setName("Um Padrao");
//			//expr.end();
//		}
//	}
//}

void Relaxacao_Linear::restricoes_demanda() {
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
					for (t = 0; t < T - TipoVigas[PackPatterns[i].tipo].tempo_cura + 1; ++t) {
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

void Relaxacao_Linear::restricoes_sequenciamento() {
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

void Relaxacao_Linear::restricoes_estoque() {
	//estoque de leftovers

	//leftovers cortados obedecem restricoes de estoque
	for (IloInt w = W; w < W + V; w++) {
		IloExpr expr(env);

		for (auto h : H_menos) {
			if (CutPatterns[h].index_barra == w)
				expr += CPLEX_y_bi[h][w];
		}

		for (auto mu : SplPatterns)
			expr += mu.tamanhos[w - W] * CPLEX_o[mu.id];

		model.add(expr <= e[w]).setName("Estoque 1");
		expr.end();
	}

	//barras novas cortadas obedecem o estoque
	for (IloInt w = 0; w < W; w++) {
		IloExpr expr1(env), expr2(env);
		for (auto h : H_menos) {
			if (CutPatterns[h].index_barra == w)
				expr1 += CPLEX_y_bi[h][w];
		}



		for (IloInt v = 0; v < V; v++) {
			for (auto h : H_mais) {
				if (CutPatterns[h].Gera_LeftOver(Gamma, v) && (CutPatterns[h].index_barra == w)) {
					expr2 += CPLEX_y_tri[h][w][v];
				}
			}
		}

		model.add(expr1 + expr2 <= e[w]).setName("Estoque 2");
		expr1.end();
		expr2.end();
	}
}

void Relaxacao_Linear::restricoes_z() {
	IloInt t, m, i;

	for (t = 0; t < T; t++)
	{
		IloExpr expr(env);

		for (m = 0; m < M; m++) {
			for (i = 1; i < P; i++)
				if (PackPatterns[i].maximal(FORMAS[m], TipoVigas[PackPatterns[i].tipo].comprimentos))
					expr += CPLEX_x[i][m][t];
			expr += CPLEX_x[0][m][t];
		}

		model.add(M*z[t] >= expr);
		expr.end();
	}
}

void Relaxacao_Linear::restricoes_continuidade() {
	IloInt m, t, i;

	for (m = 0; m < M; m++) {
		for (t = 0; t < T - 1; t++) {
			IloExpr expr1(env), expr2(env);

			for (i = 0; i < P; i++)
				if (PackPatterns[i].maximal(FORMAS[m], TipoVigas[PackPatterns[i].tipo].comprimentos) || i == 0)
					expr1 += CPLEX_x[i][m][t];
			for (i = 0; i < P; i++)
				if (PackPatterns[i].maximal(FORMAS[m], TipoVigas[PackPatterns[i].tipo].comprimentos) || i == 0)
					expr2 += CPLEX_x[i][m][t + 1];

			model.add(expr1 >= expr2);
			expr1.end();
			expr2.end();
		}
	}

}

void Relaxacao_Linear::restricoes_integracao() {

	for (IloInt gamma = 0; gamma < Gamma; gamma++)
	{
		IloExpr soma1(env), soma2(env), soma3(env), expr(env), soma4(env);

		for (IloInt w = 0; w < W; w++) {
			for (auto h : H_menos)
				if (CutPatterns[h].index_barra == w)
					soma1 += CutPatterns[h].tamanhos[gamma] * CPLEX_y_bi[h][w];
		}

		for (IloInt w = 0; w < W; w++) {
			for (IloInt v = 0; v < V; v++)
				for (auto h : H_mais)
					if (CutPatterns[h].index_barra == w && CutPatterns[h].Gera_LeftOver(Gamma, v))
						soma2 += CutPatterns[h].tamanhos[gamma] * CPLEX_y_tri[h][w][v];
		}
		for (IloInt w = W; w < W + V; w++) {
			for (auto h : H_menos)
				if (CutPatterns[h].index_barra == w)
					soma3 += CutPatterns[h].tamanhos[gamma] * CPLEX_y_bi[h][w];
		}

		for (int mu = 0; mu < O; mu++) {
			if (SplPatterns[mu].barra_gerada == gamma)
				expr += CPLEX_o[mu];

		}

		for (auto m : G[gamma]) {
			for (int t = 0; t < T; t++)
				for (int i = 1; i < P; i++)
					if (PackPatterns[i].maximal(FORMAS[m], TipoVigas[PackPatterns[i].tipo].comprimentos))
						soma4 += TipoVigas[PackPatterns[i].tipo].n_barras * CPLEX_x[i][m][t];
		}
		model.add(soma1 + soma2 + soma3 + expr == soma4).setName("Integracao");

		expr.end();
		soma1.end();
		soma2.end();
		soma3.end();
		soma4.end();
	}

}


void Relaxacao_Linear::restricao_limite() {

	IloExpr expr(env);
	IloInt estoque_leftover = 0;

	for (IloInt w = W; w < W + V; w++)
		estoque_leftover += e[w];



	for (IloInt w = 0; w < W; w++)
		for (IloInt v = 0; v < V; v++)
			for (auto h : H_mais)
				if (CutPatterns[h].Gera_LeftOver(Gamma, v))
					expr += CPLEX_y_tri[h][w][v];


	for (IloInt w = W; w < W + V; w++)
		for (auto h : H_menos)
			expr -= CPLEX_y_bi[h][w];

	for (auto mu : SplPatterns)
		for (int v = 0; v < V; v++)
			expr -= mu.tamanhos[v] * CPLEX_o[mu.id];

	U = 1000;
	model.add(expr <= U - estoque_leftover).setName("Limite de leftover");

	expr.end();
}

void Relaxacao_Linear::MontarModelo() {
	try {
		model = IloModel(env);

		iniciar_variaveis();

		CPLEX_objective_function();

		restricoes_onlyone();

		restricoes_demanda();

		restricoes_sequenciamento();

		restricoes_z();

		restricoes_continuidade();

		restricoes_estoque();

		restricoes_integracao();

		//restricao_limite();


		cplex = IloCplex(model);
		cplex.exportModel("Modelo.lp");
	}
	catch (IloException& e) {
		cerr << "Erro: " << e.getMessage() << endl;
		cout << "\nErro ilocplex" << endl;
		getchar();
		return;
	}
	catch (const exception& e) {
		cerr << "Erro: " << e.what() << endl;
		getchar();
		return;
	}
	catch (...) {
		return;
	}
}

void Relaxacao_Linear::ImprimirSolucaoArquivo() {
	ofstream resultados("resultados_rel_lin.txt", fstream::app);

	resultados << endl << nome_instancia << "," << calculo_lowerbound();
	/*resultados << "," << C << "," << M << "," << T <<
	"," << P << "," << H << "," << O << endl;*/
	try
	{
		resultados << "," << cplex.getObjValue() << "," << cplex.getValue(Makespan) << "," << cplex.getValue(Custo)
			<< "," << tempo_solucao << endl;
	}
	catch (IloException &e) {
		resultados << "," << e.getMessage() << "," << tempo_solucao << endl;
	}
	catch (...)
	{
		resultados << ",No solution," << tempo_solucao << endl;
	}

	resultados.close();
}


double Relaxacao_Linear::calculo_lowerbound()
{
	double makespan = 0.0, soma_aux = 0.0;

	//Lower bound do makespan
	for (int i = 0; i < C; i++)
	{
		for (int k = 0; k < TipoVigas[i].n_comprimentos; k++)
		{
			makespan += TipoVigas[i].tempo_cura* TipoVigas[i].comprimentos[k] * TipoVigas[i].demandas[k];
		}
	}
	for (int m = 0; m < M; m++)
	{
		soma_aux += FORMAS[m];
	}

	makespan = ceil(makespan / soma_aux);


	//Lower bound da sobra


	double melhor_sobra = INT_MAX;
	for (int gamma = 0; gamma < Gamma; gamma++)
	{
		double qtde_min_barras = 0;
		for (int i = 0; i < C; i++)
		{
			for (int k = 0; k < TipoVigas[i].n_comprimentos; k++)
			{
				qtde_min_barras += TipoVigas[i].n_barras* TipoVigas[i].comprimentos[k] * TipoVigas[i].demandas[k];
			}
		}
		qtde_min_barras = ceil(qtde_min_barras / L[gamma]);

		//min
		double melhor_padrao = INT_MAX;
		//min f_hw
		for (auto h : H_menos) {
			if (CutPatterns[h].index_barra < W)
				if (CutPatterns[h].tamanhos[gamma] > 0 && (b[CutPatterns[h].index_barra] - CutPatterns[h].cap < melhor_padrao)/ CutPatterns[h].tamanhos[gamma])
					melhor_padrao = (b[CutPatterns[h].index_barra] - CutPatterns[h].cap)/ CutPatterns[h].tamanhos[gamma];
		}
		//min alpha' f_hwv
		for (auto h : H_mais) {
			if (CutPatterns[h].tamanhos[gamma] > 0 && Parameter_alpha_1*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap < melhor_padrao)/ CutPatterns[h].tamanhos[gamma])
				melhor_padrao = Parameter_alpha_1*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap)/ CutPatterns[h].tamanhos[gamma];
		}

		//min alpha'' fhw
		for (auto h : H_menos) {
			if (CutPatterns[h].index_barra >= W)
				if (CutPatterns[h].tamanhos[gamma] > 0 && Parameter_alpha_2*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap < melhor_padrao)/ CutPatterns[h].tamanhos[gamma])
					melhor_padrao = Parameter_alpha_2*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap)/ CutPatterns[h].tamanhos[gamma];
		}

		for (auto mu : SplPatterns) {
			if (mu.barra_gerada == gamma && Parameter_alpha_2 * mu.folga < melhor_padrao)
				melhor_padrao = Parameter_alpha_2 * mu.folga;
		}


		if (melhor_sobra > qtde_min_barras*melhor_padrao)
			melhor_sobra = qtde_min_barras*melhor_padrao;
	}


	return makespan + melhor_sobra;
}