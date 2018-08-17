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
	Tipos_de_Viga.resize(C);
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

		instancia >> Tipos_de_Viga[i].tempo_cura
			>> Tipos_de_Viga[i].n_comprimentos
			>> Tipos_de_Viga[i].n_barras;
		Tipos_de_Viga[i].demandas.resize(Tipos_de_Viga[i].n_comprimentos);
		Tipos_de_Viga[i].comprimentos.resize(Tipos_de_Viga[i].n_comprimentos);

		if (Tipos_de_Viga[i].n_comprimentos > Maior_Qc)
			Maior_Qc = Tipos_de_Viga[i].n_comprimentos;



		//Ler comprimentos
		for (int k = 0; k < Tipos_de_Viga[i].n_comprimentos; k++) {
			instancia >> Tipos_de_Viga[i].comprimentos[k];
			if (Tipos_de_Viga[i].comprimentos[k] < Menor_tamanho[i])
				Menor_tamanho[i] = Tipos_de_Viga[i].comprimentos[k];
		}

		//Demandas dos comprimentos
		for (int k = 0; k < Tipos_de_Viga[i].n_comprimentos; k++)
			instancia >> Tipos_de_Viga[i].demandas[k];

	}

	b.resize(W + V); //Alocar vetor com os tamanhos únicos de barras

	instancia.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de corte*/
	ler_pat_cut >> H;
	Padroes_de_Corte.resize(H);
	for (int h = 0; h < H; h++){
		double soma_cap = 0;

		Padroes_de_Corte[h].index_pat = h;
		ler_pat_cut >> Padroes_de_Corte[h].index_barra;
		Padroes_de_Corte[h].tamanhos.resize(Gamma + V);
		for (int i = 0; i < Gamma + V; i++){
			ler_pat_cut >> Padroes_de_Corte[h].tamanhos[i];
			if (i < Gamma)
				soma_cap += Padroes_de_Corte[h].tamanhos[i] * L[i];
			else
				soma_cap += Padroes_de_Corte[h].tamanhos[i] * b[i];
		}
		Padroes_de_Corte[h].cap = soma_cap;
	}
	ler_pat_cut.close();

	/*Leitura dos dados e alocação*/
	/*Padroes de empacotamento*/
	ler_pat_pack >> P;
	Padroes_de_Empacotamento.resize(P);
	for (int i = 0; i < P; i++) {
		Padroes_de_Empacotamento[i].id = i;
		ler_pat_cut >> Padroes_de_Empacotamento[i].tipo;

		Padroes_de_Empacotamento[i].n_comprimentos
			= Tipos_de_Viga[Padroes_de_Empacotamento[i].tipo].n_comprimentos;


		Padroes_de_Empacotamento[i].tamanhos.resize(Padroes_de_Empacotamento[i].n_comprimentos);

		int contar_cobertos = 0;
		double soma = 0;
		for (int k = 0; k < Padroes_de_Empacotamento[i].n_comprimentos; k++) {
			ler_pat_cut >> Padroes_de_Empacotamento[i].tamanhos[k];
			soma += Tipos_de_Viga[Padroes_de_Empacotamento[i].tipo].comprimentos[k]
				* Padroes_de_Empacotamento[i].tamanhos[k];

			if (Padroes_de_Empacotamento[i].tamanhos[k] > 0)
				contar_cobertos++;
		}
		Padroes_de_Empacotamento[i].n_cobertos = contar_cobertos;
		Padroes_de_Empacotamento[i].cap = soma;

	}
	ler_pat_pack.close();


	Maior_Valor_nos_Padroes = 0; //Será usado para calcular os upper bounds das variáveis de decisão

	for (int h = 0; h < H; h++)
	{
		for (int i = 0; i < Gamma + V; i++) {
			if (Padroes_de_Corte[h].tamanhos[i] > Maior_Valor_nos_Padroes)
				Maior_Valor_nos_Padroes = Padroes_de_Corte[h].tamanhos[i];
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
