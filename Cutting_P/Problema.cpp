#include "Problema.h"



Problema::Problema(const char* filename) {
	//Criar nome dos arquivos de padroes com base no nome da instância
	Parameter_alpha_1 = 2;
	Parameter_alpha_2 = 0;

	string padroes_cut, padroes_pack, padroes_spl;
	nome_instancia = filename;
	stringstream ss,ss2,ss3;
	ss << filename << ".cutpat";
	padroes_cut = ss.str();
	ss2 << filename << ".pat";
	padroes_pack = ss2.str();
	ss3 << filename << ".spl";
	padroes_spl = ss3.str();

	ifstream instancia(filename, ifstream::in);
	ifstream ler_pat_pack(padroes_pack, ifstream::in);
	ifstream ler_pat_cut(padroes_cut, ifstream::in);
	ifstream ler_pat_spl(padroes_spl, ifstream::in);

	//Testar se as instâncias foram encontradas
	bool fail = false;
	if (instancia.fail()) {
		cerr << "     Arquivo \"" << filename << "\" nao encontrado." << endl;
		fail = true;
	}
	if (ler_pat_cut.fail())
	{
		cerr << "     Arquivo \"" << padroes_cut << "\" nao encontrado." << endl;
		fail = true;
	}
	if (ler_pat_pack.fail())
	{
		cerr << "     Arquivo \"" << padroes_pack << "\" nao encontrado." << endl;
		fail = true;
	}
	if (ler_pat_spl.fail())
	{
		cerr << "     Arquivo \"" << padroes_pack << "\" nao encontrado." << endl;
		fail = true;
	}
	if (fail) {
		system("pause");
		exit(0);
	}

	/*Leitura dos dados e alocação*/
	/*Problema em si*/

	
	instancia >> C >> M >> T; //Só usa para padroes de empacotamento
	cout << C << " " << M << " " << T << endl;
	instancia >> W >> V; //Só Usa para padroes de corte
	cout << W << " " << V << endl;

	L.resize(M); //vai ser vetor unique
	FORMAS.resize(M);
	TipoVigas.resize(C);
	Maior_Forma = 0;
	Menor_Forma = INT_MAX;

	//Ler vetor das formas e seus tamanhos
	for (int i = 0; i < M; i++) {
		instancia >> FORMAS[i];
		cout << FORMAS[i] << " ";
		if (FORMAS[i] > Maior_Forma)
			Maior_Forma = FORMAS[i];
		if ((FORMAS[i] < Menor_Forma))
			Menor_Forma = FORMAS[i];
	}
	cout << endl << endl;
	L = FORMAS;
	sort(L.begin(), L.end());
	/*Calculo do numero de tamanhos de forma diferente e seu vetor*/
	auto iterador_auxiliar = unique(L.begin(), L.end());//Tirar duplicados
	L.resize(distance(L.begin(), iterador_auxiliar));//Mudar o tamanho
	Gamma = L.size();//Gamma é o tamanho do vetor de comprimentos únicos de forma

					 /*Para cada tipo de viga ler os n_tamanhos, tamanhos, demandas, cura e n_barras
					 necessárias*/
	Maior_Qc = 0;
	Menor_tamanho.resize(C);
	for (int i = 0; i < C; i++) {
		Menor_tamanho[i] = 100;

		instancia >> TipoVigas[i].tempo_cura
			>> TipoVigas[i].n_comprimentos
			>> TipoVigas[i].n_barras;

		cout << endl << endl << TipoVigas[i].tempo_cura
			<< " " << TipoVigas[i].n_comprimentos
			<< " " << TipoVigas[i].n_barras << endl;

		TipoVigas[i].demandas.resize(TipoVigas[i].n_comprimentos);
		TipoVigas[i].comprimentos.resize(TipoVigas[i].n_comprimentos);

		if (TipoVigas[i].n_comprimentos > Maior_Qc)
			Maior_Qc = TipoVigas[i].n_comprimentos;


		//Ler comprimentos
		for (int k = 0; k < TipoVigas[i].n_comprimentos; k++) {
			instancia >> TipoVigas[i].comprimentos[k];
			cout << TipoVigas[i].comprimentos[k] << " ";
			if (TipoVigas[i].comprimentos[k] < Menor_tamanho[i])
				Menor_tamanho[i] = TipoVigas[i].comprimentos[k];
		}

		//Demandas dos comprimentos
		for (int k = 0; k < TipoVigas[i].n_comprimentos; k++)
			instancia >> TipoVigas[i].demandas[k];

	}

	b.resize(W + V); //Alocar vetor com os tamanhos únicos de barras
	e.resize(W + V);

	cout << endl << endl;
	//Ler tamanhos
	Maior_Barra = 0;
	for (int i = 0; i < W + V; i++) {
		instancia >> b[i];
		cout << b[i] << " ";
		if (b[i] > Maior_Barra)
			Maior_Barra = b[i];
	}
	//Ler estoque
	cout << endl;
	for (int i = 0; i < W + V; i++) {
		instancia >> e[i];
		cout << e[i] << " ";
	}
	instancia >> SPL_epsilon;
	cout << endl << SPL_epsilon;
	instancia.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de traspasse*/

	ler_pat_spl >> O;
	SplPatterns.resize(O);
	for (int mu = 0; mu < O; mu++){
		SplPatterns[mu].id = mu;
		ler_pat_spl >> SplPatterns[mu].barra_gerada >>
			SplPatterns[mu].folga;
		SplPatterns[mu].tamanhos.resize(V);
		for (int v = 0; v < V; v++)
			ler_pat_spl >> SplPatterns[mu].tamanhos[v];
	}
	ler_pat_spl.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de corte*/
	ler_pat_cut >> H;
	CutPatterns.resize(H);
	for (int h = 0; h < H; h++){

		CutPatterns[h].index_pat = h;
		ler_pat_cut >> CutPatterns[h].index_barra;
		ler_pat_cut >> CutPatterns[h].cap;

		CutPatterns[h].tamanhos.resize(Gamma + V);
		for (int i = 0; i < Gamma + V; i++) 
			ler_pat_cut >> CutPatterns[h].tamanhos[i];

		
		
	}
	ler_pat_cut.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de empacotamento*/
	ler_pat_pack >> P;
	PackPatterns.resize(P);
	for (int i = 0; i < P; i++) {
		PackPatterns[i].id = i;
		ler_pat_pack >> PackPatterns[i].tipo;
		ler_pat_pack >> PackPatterns[i].cap;
		PackPatterns[i].n_comprimentos
			= TipoVigas[PackPatterns[i].tipo].n_comprimentos;


		PackPatterns[i].tamanhos.resize(PackPatterns[i].n_comprimentos);

		int contar_cobertos = 0;

		for (int k = 0; k < PackPatterns[i].n_comprimentos; k++) {
			ler_pat_pack >> PackPatterns[i].tamanhos[k];
			if (PackPatterns[i].tamanhos[k] > 0)
				contar_cobertos++;
		}
		PackPatterns[i].n_cobertos = contar_cobertos;

	}
	ler_pat_pack.close();


	Maior_Valor_nos_Padroes = 0; //Será usado para calcular os upper bounds das variáveis de decisão

	for (int h = 0; h < H; h++)
	{
		for (int i = 0; i < Gamma + V; i++) {
			if (CutPatterns[h].tamanhos[i] > Maior_Valor_nos_Padroes)
				Maior_Valor_nos_Padroes = CutPatterns[h].tamanhos[i];
		}
	}
	
	cout << endl;
	instancia.close();

	for (int i = 0; i < 50; i++) cout << "_";
	cout << endl;
	cout << "\t \t Leitura da instancia completa" << endl;
	for (int i = 0; i < 50; i++) cout << "_";
	cout << endl;


}


Problema::~Problema()
{
}
