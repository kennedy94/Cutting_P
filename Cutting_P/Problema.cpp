#include "Problema.h"



Problema::Problema(const char* filename) {
	//Criar nome dos arquivos de padroes com base no nome da instância
	string padroes_cut, padroes_pack;
	nome_instancia = filename;
	stringstream ss;
	ss << filename << ".cutpat";
	padroes_cut = ss.str();
	ss.clear();
	ss << filename << ".pat";
	padroes_pack = ss.str();

	ifstream instancia(filename, ifstream::in);
	ifstream ler_pat_pack(padroes_pack, ifstream::in);
	ifstream ler_pat_cut(padroes_cut, ifstream::in);

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
	if (fail) {
		system("pause");
		exit(0);
	}

	/*Leitura dos dados e alocação*/
	/*Problema em si*/

	
	instancia >> C >> M >> T; //Só usa para padroes de empacotamento
	instancia >> W >> V; //Só Usa para padroes de corte

	L.resize(M); //vai ser vetor unique
	FORMAS.resize(M);
	TipoVigas.resize(C);
	Maior_Forma = 0;
	Menor_Forma = INT_MAX;

	//Ler vetor das formas e seus tamanhos
	for (int i = 0; i < M; i++) {
		instancia >> FORMAS[i];
		if (FORMAS[i] > Maior_Forma)
			Maior_Forma = FORMAS[i];
		if ((FORMAS[i] < Menor_Forma))
			Menor_Forma = FORMAS[i];
	}

	L = FORMAS;
	/*Calculo do numero de tamanhos de forma diferente e seu vetor*/
	auto iterador_auxiliar = unique(L.begin(), L.end());//Tirar duplicados
	FORMAS.resize(distance(L.begin(), iterador_auxiliar));//Mudar o tamanho
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
		TipoVigas[i].demandas.resize(TipoVigas[i].n_comprimentos);
		TipoVigas[i].comprimentos.resize(TipoVigas[i].n_comprimentos);

		if (TipoVigas[i].n_comprimentos > Maior_Qc)
			Maior_Qc = TipoVigas[i].n_comprimentos;



		//Ler comprimentos
		for (int k = 0; k < TipoVigas[i].n_comprimentos; k++) {
			instancia >> TipoVigas[i].comprimentos[k];
			if (TipoVigas[i].comprimentos[k] < Menor_tamanho[i])
				Menor_tamanho[i] = TipoVigas[i].comprimentos[k];
		}

		//Demandas dos comprimentos
		for (int k = 0; k < TipoVigas[i].n_comprimentos; k++)
			instancia >> TipoVigas[i].demandas[k];

	}

	b.resize(W + V); //Alocar vetor com os tamanhos únicos de barras
	e.resize(W + V);

	//Ler tamanhos
	for (int i = 0; i < W + V; i++)
		instancia >> b[i];
	//Ler estoque
	for (int i = 0; i < W + V; i++)
		instancia >> e[i];

	instancia.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de corte*/
	ler_pat_cut >> H;
	CutPatterns.resize(H);
	for (int h = 0; h < H; h++){
		double soma_cap = 0;

		CutPatterns[h].index_pat = h;
		ler_pat_cut >> CutPatterns[h].index_barra;
		CutPatterns[h].tamanhos.resize(Gamma + V);
		for (int i = 0; i < Gamma + V; i++){
			ler_pat_cut >> CutPatterns[h].tamanhos[i];
			if (i < Gamma)
				soma_cap += CutPatterns[h].tamanhos[i] * L[i];
			else
				soma_cap += CutPatterns[h].tamanhos[i] * b[i];
		}
		CutPatterns[h].cap = soma_cap;
	}
	ler_pat_cut.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de empacotamento*/
	ler_pat_pack >> P;
	PackPatterns.resize(P);
	for (int i = 0; i < P; i++) {
		PackPatterns[i].id = i;
		ler_pat_cut >> PackPatterns[i].tipo;

		PackPatterns[i].n_comprimentos
			= TipoVigas[PackPatterns[i].tipo].n_comprimentos;


		PackPatterns[i].tamanhos.resize(PackPatterns[i].n_comprimentos);

		int contar_cobertos = 0;
		double soma = 0;
		for (int k = 0; k < PackPatterns[i].n_comprimentos; k++) {
			ler_pat_cut >> PackPatterns[i].tamanhos[k];
			soma += TipoVigas[PackPatterns[i].tipo].comprimentos[k]
				* PackPatterns[i].tamanhos[k];

			if (PackPatterns[i].tamanhos[k] > 0)
				contar_cobertos++;
		}
		PackPatterns[i].n_cobertos = contar_cobertos;
		PackPatterns[i].cap = soma;

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


	for (int i = 0; i < 10; i++) cout << "_";
	cout << endl;
	cout << "\t \t Leitura da instancia completa" << endl;
	for (int i = 0; i < 10; i++) cout << "_";
	cout << endl;


}


Problema::~Problema()
{
}
