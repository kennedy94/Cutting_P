#include "Modelo_Cplex.h"


bool GreaterThanZero(int i) { return (i > 0); }

void Modelo_Cplex::iniciar_variaveis(){

	char strnum[30];

	z = IloBoolVarArray(env, T);
	for (int t = 0; t < T; t++) {
		sprintf(strnum, "z(%d)", t);
		z[t].setName(strnum);
	}


	CPLEX_o = IloIntVarArray(env, O, 0, 10000000);

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

	CPLEX_y_bi = IloArray<IloIntVarArray>(env, H);
	for (int h = 0; h < H; h++) {
		CPLEX_y_bi[h] = IloIntVarArray(env, W + V, 0, 10000000);
		for (IloInt w = 0; w < W+V; w++)
		{
			sprintf(strnum, "y(%d,%d)", h, w);
			CPLEX_y_bi[h][w].setName(strnum);
		}
	}
	
	
	
	CPLEX_y_tri = IloArray<IloArray<IloIntVarArray>>(env, H);
	for (int h = 0; h < H; h++) {
		CPLEX_y_tri[h] = IloArray<IloIntVarArray>(env, W);
		for (IloInt w = 0; w < W; w++) {
			CPLEX_y_tri[h][w] = IloIntVarArray(env, V, 0, 10000000);
			for (IloInt v = 0; v < V; v++)
			{
				sprintf(strnum, "y(%d,%d,%d)", h,w,v);
				CPLEX_y_tri[h][w][v].setName(strnum);
			}
		}
	}
	
	CPLEX_x = IloArray<IloArray<IloBoolVarArray>>(env, P);
	for (IloInt i = 0; i < P; i++) {
		CPLEX_x[i] = IloArray<IloBoolVarArray>(env, M);
		for (IloInt m = 0; m < M; m++) {
			CPLEX_x[i][m] = IloBoolVarArray(env, T);
			for (IloInt t = 0; t < T; t++)
			{
				sprintf(strnum, "x(%d,%d,%d)", i, m, t);
				CPLEX_x[i][m][t].setName(strnum);
			}
		}
	}
}

void Modelo_Cplex::resolver_inteira() {
	cplex = IloCplex(model);

	try
	{
		//cplex.setParam(IloCplex::PreInd, 0);
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

Modelo_Cplex::~Modelo_Cplex()
{
}

void Modelo_Cplex::CPLEX_objective_function(){
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



//SOS
//void Modelo_Cplex::restricoes_onlyone() {
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

		for (auto h : H_menos) {
			if(CutPatterns[h].index_barra == w)
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

void Modelo_Cplex::restricoes_z() {
	IloInt t, m, i;

	for (t = 0; t < T; t++)
	{
		IloExpr expr(env);

		for (m = 0; m < M; m++) {
			for (i = 1; i < P; i++)
				if(PackPatterns[i].maximal(FORMAS[m],TipoVigas[PackPatterns[i].tipo].comprimentos))
					expr += CPLEX_x[i][m][t];
			expr += CPLEX_x[0][m][t];
		}

		model.add(M*z[t] >= expr);
		expr.end();
	}
}

void Modelo_Cplex::restricoes_continuidade() {
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

void Modelo_Cplex::restricoes_integracao() {

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
			if(SplPatterns[mu].barra_gerada == gamma)
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


void Modelo_Cplex::restricao_limite() {

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

void Modelo_Cplex::MontarModelo() {
	try {
		model = IloModel(env);

		iniciar_variaveis();

		//CPLEX_objective_function_transformada();

		CPLEX_objective_function_transformada();

		restricoes_onlyone();

		restricoes_demanda();

		restricoes_sequenciamento();

		restricoes_z();

		restricoes_continuidade();

		restricoes_estoque();

		restricoes_integracao();

		//restricao_limite();

		
		cplex = IloCplex(model);
		//cplex.exportModel("Modelo.lp");
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

void Modelo_Cplex::ImprimirSolucaoArquivo(){
	ofstream resultados("resultados.txt", fstream::app);

	resultados << endl << nome_instancia;
	/*resultados << "," << C << "," << M << "," << T << 
		"," << P << "," << H << "," << O << endl;*/
	try
	{
		resultados << "," << cplex.getObjValue() << "," << cplex.getValue(Makespan) << "," << cplex.getValue(Custo) << "," <<
			cplex.getNnodes() << "," << cplex.getMIPRelativeGap() << "," << tempo_solucao << endl;
	}
	catch (...)
	{
		resultados << ",No solution," << tempo_solucao <<  endl;
	}
	
	resultados.close();
}

void Modelo_Cplex::ImprimirSolucao() {
	int contador = 1;
	for (int t = 0; t < T; t++)
		for (int m = 0; m < M; m++) 
			for (int i = 0; i < P; i++) 
				if (PackPatterns[i].maximal(FORMAS[m], TipoVigas[PackPatterns[i].tipo].comprimentos) || i == 0)
					if (cplex.isExtracted(CPLEX_x[i][m][t]) && cplex.getValue(CPLEX_x[i][m][t]) > 0.00001) {
						cout << "     x(" << i << "," << m << "," << t << ") = " << cplex.getValue(CPLEX_x[i][m][t]) << "		";
						if (contador % 2 == 0)
							cout << endl;
						contador++;
					}

	for (int h = 0; h < H; h++) 
		for (int w = 0; w < W; w++) 
			for (int v = 0; v < V; v++) 
				if (cplex.isExtracted(CPLEX_y_tri[h][w][v]) && cplex.getValue(CPLEX_y_tri[h][w][v]) > 0) {
					cout << "     y(" << h << "," << w << "," << v << ") = " << cplex.getValue(CPLEX_y_tri[h][w][v]) << "		";
					if (contador % 2 == 0)
						cout << endl;
					contador++;
				}

	for (int h = 0; h < H; h++) 
		for (int w = 0; w < W; w++) 
			if (cplex.isExtracted(CPLEX_y_bi[h][w]) && cplex.getValue(CPLEX_y_bi[h][w]) > 0) {
				cout << "     y(" << h << "," << w << ") = " << cplex.getValue(CPLEX_y_bi[h][w]) << "		";
				if (contador % 2 == 0)
					cout << endl;
				contador++;
			}
	

	for (auto h : SplPatterns)
		if (cplex.isExtracted(CPLEX_o[h.id]) && cplex.getValue(CPLEX_o[h.id]) > 0)
			cout << "o(" << h.id << ") = " << cplex.getValue(CPLEX_o[h.id]) << endl;

	cout << endl << "Objetivo -> " << cplex.getObjValue() << endl;
	
}


void Modelo_Cplex::ImprimirGantt() {
	char xu[100];
	strcpy(xu, nome_instancia);
	strcat(xu, ".solu");
	ofstream txtsolu;
	txtsolu.open(xu, fstream::trunc);

	txtsolu << endl;

	txtsolu << 1 << "," << 0 << "," << T << ",Type 0" << endl;
	for (int m = 0; m < M; m++) {
		bool usou = false;
		for (int t = 0; t < T; t++)
			for (int i = 1; i < P; i++)
				if (PackPatterns[i].maximal(FORMAS[m], TipoVigas[PackPatterns[i].tipo].comprimentos)
					&& cplex.isExtracted(CPLEX_x[i][m][t])
					&& cplex.getValue(CPLEX_x[i][m][t]) == 1) {
					txtsolu << m + 1 << "," << t + 0.01 << "," << t + TipoVigas[PackPatterns[i].tipo].tempo_cura - 0.01 << ",Type " << PackPatterns[i].tipo + 1 << endl;
					usou = true;
				}
		if (!usou)
			txtsolu << m + 1 << "," << 0 << "," << T << ",Type 0" << endl;
	}
	txtsolu.close();

	string comando_system;

	stringstream ss;
	ss << " python imprimir_gantt.py " << xu;
	comando_system = ss.str();
	system(comando_system.c_str());
	cout << "Grafico de Gantt Gerado! :)" << endl;
}


void Modelo_Cplex::PlotarBarras() {
	char xu[100];
	strcpy(xu, nome_instancia);
	strcat(xu, ".barras");
	ofstream txtsolu;
	txtsolu.open(xu, fstream::trunc);

	vector<int> UsedCutPatterns(H);
	for (int h = 0; h < H; h++)
		UsedCutPatterns[h] = 0;


	for (int h = 0; h < H; h++)
		for (int w = 0; w < W; w++)
			for (int v = 0; v < V; v++)
				if (cplex.isExtracted(CPLEX_y_tri[h][w][v]) && cplex.getValue(CPLEX_y_tri[h][w][v]) > 0)
					UsedCutPatterns[h] = cplex.getValue(CPLEX_y_tri[h][w][v]);

	for (int h = 0; h < H; h++)
		for (int w = 0; w < W; w++)
			if (cplex.isExtracted(CPLEX_y_bi[h][w]) && cplex.getValue(CPLEX_y_bi[h][w]) > 0)
				UsedCutPatterns[h] = cplex.getValue(CPLEX_y_bi[h][w]);


	txtsolu << 1 << "," << 0 << "," << Maior_Barra << ",Type 0" << endl;

	int contador = 0;
	for (int h = 0; h < H; h++){
		if (UsedCutPatterns[h] > 0) {

			double aux = 0;
			contador++;
			

			txtsolu << contador << "," << CutPatterns[h].cap + 0.1 << "," << b[CutPatterns[h].index_barra] + 0.1 << ",Type 3" << endl;


			for (int w = 0; w < Gamma+V; w++){
				for (int i = 0; i < CutPatterns[h].tamanhos[w]; i++){
					if (w < Gamma) {
						txtsolu << contador << "," << aux + 0.1 << "," << aux + L[w] - 0.1 << ",Type 1" << endl;
						aux += L[w];
					}
					else {
						txtsolu << contador << "," << aux + 0.1 << "," << aux + b[W + w - Gamma] - 0.1 << ",Type 2" << endl;
						aux += b[W + w - Gamma];
					}
				}
			}
		}
	}

	for (auto mu : SplPatterns)
	{
		if (cplex.isExtracted(CPLEX_o[mu.id]) && cplex.getValue(CPLEX_o[mu.id]) > 0) {
			double aux = 0;
			contador++;
			double soma = 0;
			
			for (int i = 0; i < V; i++)
				soma += mu.tamanhos[i]* b[W + i];

			//txtsolu << contador << "," << 0
			//	<< "," << soma - SPL_epsilon<< ",Type 4" << endl;

			txtsolu << contador << "," <<  L[mu.barra_gerada]
				<< "," << L[mu.barra_gerada] + mu.folga - SPL_epsilon << ",Type 3" << endl;

			for (auto i : mu.tamanhos) {
				for (int v = 0; v < i; v++) {
					if(aux == 0)
						txtsolu << contador << "," << 0 <<
						"," << b[W + v] - SPL_epsilon / 2 << ",Type 4" << endl;
					if (aux > 0) {
						txtsolu << contador << "," << aux << ","
							<< aux + SPL_epsilon << ",Emenda" << endl;
						txtsolu << contador << "," << aux + SPL_epsilon <<
							"," << L[mu.barra_gerada] << ",Type 4" << endl;


					}
					aux += b[W + v] - SPL_epsilon / 2;
				}
			}

		}
	}


	txtsolu.close();

	string comando_system;

	stringstream ss;
	ss << " python gerarbarras.py " << xu;
	comando_system = ss.str();
	system(comando_system.c_str());
	cout << "Barras Gerado! :)" << endl;

}



void Modelo_Cplex::CPLEX_objective_function_transformada() {
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
	Makespan = IloNumVar(env, 0, T);
	Custo = IloNumVar(env, 0, 99999999);

	model.add(costSum == Makespan);
	model.add(Custo == soma1 + Parameter_alpha_1*soma2 + Parameter_alpha_2*(soma3 + expr));
	model.add(IloMinimize(env, Makespan/LB_makespan() + Custo/LB_custo()));

	soma1.end();
	soma2.end();
	soma3.end();
	expr.end();
	costSum.end();
}


IloNum Modelo_Cplex::LB_makespan()
{
	IloNum makespan = 0.0, soma_aux = 0.0;

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

	return makespan;
}

IloNum Modelo_Cplex::LB_custo()
{
	IloNum melhor_sobra = (double)INT_MAX;
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
		double melhor_padrao = (double)INT_MAX;
		//min f_hw
		for (auto h : H_menos) {
			if (CutPatterns[h].index_barra < W)
				if (CutPatterns[h].tamanhos[gamma] > 0 && (b[CutPatterns[h].index_barra] - CutPatterns[h].cap < melhor_padrao) / CutPatterns[h].tamanhos[gamma])
					melhor_padrao = (b[CutPatterns[h].index_barra] - CutPatterns[h].cap) / CutPatterns[h].tamanhos[gamma];
		}
		//min alpha' f_hwv
		for (auto h : H_mais) {
			if (CutPatterns[h].tamanhos[gamma] > 0 && Parameter_alpha_1*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap < melhor_padrao) / CutPatterns[h].tamanhos[gamma])
				melhor_padrao = Parameter_alpha_1*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap) / CutPatterns[h].tamanhos[gamma];
		}

		//min alpha'' fhw
		for (auto h : H_menos) {
			if (CutPatterns[h].index_barra >= W)
				if (CutPatterns[h].tamanhos[gamma] > 0 && Parameter_alpha_2*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap < melhor_padrao) / CutPatterns[h].tamanhos[gamma])
					melhor_padrao = Parameter_alpha_2*(b[CutPatterns[h].index_barra] - CutPatterns[h].cap) / CutPatterns[h].tamanhos[gamma];
		}

		for (auto mu : SplPatterns) {
			if (mu.barra_gerada == gamma && Parameter_alpha_2 * mu.folga < melhor_padrao)
				melhor_padrao = Parameter_alpha_2 * mu.folga;
		}


		if (melhor_sobra > qtde_min_barras*melhor_padrao)
			melhor_sobra = qtde_min_barras*melhor_padrao;
	}


	return melhor_sobra;

}