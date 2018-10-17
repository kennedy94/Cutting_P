#include "Heuristica.h"



struct Utilizacao {
	int idx, util;
	bool operator<(const Utilizacao& rhs) const {
		return util < rhs.util;
	}
};





double Heuristica::fitness(individuo solu)
{
	vector<int> UtilizacaoFormas(M, 0);
	//vector<vector<int>> PadroesAssociados(M);

	int Makespan;
	double 
		Sobra1 = 0.0,
		Sobra2 = 0.0,
		Sobra3 = 0.0;
	

	
	vector<Utilizacao> UtilFormas(M);
	for (int m = 0; m < M; m++) {
		UtilFormas[m].idx = m;
		UtilFormas[m].util = 0;
	}


	for (int i = 0; i < P-1; i++)
	{
		if (solu.n_vezes[i] == 0)
			continue;

		int tipo_da_forma = -1;
		for (int gamma = 0; gamma < Gamma; gamma++) {
			if (PackPatterns[solu.ind[i]].maximal(L[gamma], TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos)) {
				tipo_da_forma = gamma;
				break;
			}
		}

		for (int vezes = 0; vezes < solu.n_vezes[i]; vezes++) {
			sort(UtilFormas.begin(), UtilFormas.end());
			//forma atual recebe argmin da utilização que acomoda o padrão
			int k_esimo = 0;
			int FormaAtual, k_escolhido;
			do
			{	
				FormaAtual = UtilFormas[k_esimo].idx;
				k_escolhido = k_esimo;
				k_esimo++;
			} while (FORMAS[FormaAtual] != L[tipo_da_forma]);

			for (int m = 0; m < M; m++) 
				if (UtilFormas[m].idx == FormaAtual) 
					UtilFormas[k_escolhido].util += TipoVigas[PackPatterns[solu.ind[i]].tipo].tempo_cura;
		}
	}
	Utilizacao valor = *max_element(UtilFormas.begin(), UtilFormas.end());
	Makespan = valor.util;

	for (int i = P - 1; i < P - 1 + H; i++){
		if (CutPatterns[solu.ind[i]].index_barra < W & !CutPatterns[solu.ind[i]].Gera_LeftOvers(Gamma, V))
			Sobra1 += solu.n_vezes[i] * (b[CutPatterns[solu.ind[i]].index_barra] - CutPatterns[solu.ind[i]].cap);

		if (CutPatterns[solu.ind[i]].index_barra < W & CutPatterns[solu.ind[i]].Gera_LeftOvers(Gamma, V))
			Sobra2 += solu.n_vezes[i] * (b[CutPatterns[solu.ind[i]].index_barra] - CutPatterns[solu.ind[i]].cap);

		if(CutPatterns[solu.ind[i]].index_barra >= W)
			Sobra3 += solu.n_vezes[i] * (b[CutPatterns[solu.ind[i]].index_barra] - CutPatterns[solu.ind[i]].cap);
	}
	for (int i = P - 1 + H; i < P - 1 + H + O; i++)
		Sobra3 += solu.n_vezes[i] * (SplPatterns[solu.ind[i]].folga);
	
	return Makespan + Sobra1 + Parameter_alpha_1*Sobra2 + Parameter_alpha_2*Sobra3;
}









list<individuo> Heuristica::cruzar(individuo pai, individuo mae)
{
	list<individuo> filhos;

	/*default_random_engine generator;
	generator.seed(time(NULL));
	uniform_int_distribution<int> bloco1(0, P - 2);
	uniform_int_distribution<int> bloco2(P - 1, P - 2 + H);

	int escollhidoPack = bloco1(generator),
		escollhidoCut = bloco2(generator);*/


	//----------------------GERAÇÃO DO FILHO1
	individuo filho = pai;

	for (int i = 0; i < P - 1 + H + O; i++) {
		filho.ind[i] = -1;
		filho.n_vezes[i] = -1;
	}


	//CROSSOVER INDICES DE PADRÕES DE EMPACOTAMENTO
	for (int i = 0; i < floor((double)(P - 1)/2); i++) {
		filho.ind[i] = pai.ind[i];
	}

	vector<int> nao_escolhidos;
	for (int i = floor((double)(P - 1) / 2); i < P - 1; i++)
		nao_escolhidos.push_back(pai.ind[i]);

	for (int i = floor((double)(P - 1) / 2); i < P - 1; i++) {
		if (find(nao_escolhidos.begin(), nao_escolhidos.end(), mae.ind[i]) != nao_escolhidos.end()) {
			filho.ind[i] = mae.ind[i];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
		else {
			default_random_engine generator;
			generator.seed(time(NULL));
			uniform_int_distribution<int> distribution(0, nao_escolhidos.size()-1);

			int escolhido = distribution(generator);
			
			filho.ind[i] = nao_escolhidos[escolhido];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
	}

	//CROSSOVER INDICES DE PADRÕES DE CORTE
	for (int i = P - 1; i < P - 1 + floor((double)H/2); i++) {
		filho.ind[i] = pai.ind[i];
	}


	for (int i = P - 1 + floor((double)H / 2); i < P - 1 + H; i++)
		nao_escolhidos.push_back(pai.ind[i]);

	for (int i = P - 1 + floor((double)H / 2); i < P - 1 + H; i++) {
		if (find(nao_escolhidos.begin(), nao_escolhidos.end(), mae.ind[i]) != nao_escolhidos.end()) {
			filho.ind[i] = mae.ind[i];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
		else {
			default_random_engine generator;
			generator.seed(time(NULL));
			uniform_int_distribution<int> distribution(0, nao_escolhidos.size() - 1);

			int escolhido = distribution(generator);

			filho.ind[i] = nao_escolhidos[escolhido];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
	}

	//CROSSOVER INDICES DE PADRÕES DE TRASPASSE
	for (int i = P - 1 + H; i <P - 1 + H + floor((double)O/2); i++) {
		filho.ind[i] = pai.ind[i];
	}
	for (int i = P - 1 + H + floor((double)O / 2); i < P - 1 + H + O; i++)
		nao_escolhidos.push_back(pai.ind[i]);

	for (int i = P - 1 + H + floor((double)O / 2); i < P - 1 + H + O; i++) {
		if (find(nao_escolhidos.begin(), nao_escolhidos.end(), mae.ind[i]) != nao_escolhidos.end()) {
			filho.ind[i] = mae.ind[i];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
		else {
			default_random_engine generator;
			generator.seed(time(NULL));
			uniform_int_distribution<int> distribution(0, nao_escolhidos.size() - 1);

			int escolhido = distribution(generator);

			filho.ind[i] = nao_escolhidos[escolhido];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
	}

	//CALCULAR NOVAS QUANTIDADES PELAS MÉDIAS
	for (int i = 0; i < P - 1; i++) {
		int	valor1 = -1,
			valor2 = -1;
		for (int j = 0; j < P - 1; j++) {
			if (pai.ind[j] == filho.ind[i])
				valor1 = pai.n_vezes[j];
			if (mae.ind[j] == filho.ind[i])
				valor2 = mae.n_vezes[j];
			if (valor1 != -1 && valor2 != -1)
				break;
		}
		filho.n_vezes[i] = ceil((double)(valor1 + valor2) / 2);
	}

	for (int i = P - 1; i < P - 1 + H; i++) {
		int	valor1 = -1,
			valor2 = -1;
		for (int j = P - 1; j < P - 1 + H; j++) {
			if (pai.ind[j] == filho.ind[i])
				valor1 = pai.n_vezes[j];
			if (mae.ind[j] == filho.ind[i])
				valor2 = mae.n_vezes[j];
			if (valor1 != -1 && valor2 != -1)
				break;
		}
		filho.n_vezes[i] = ceil((double)(valor1 + valor2) / 2);
	}

	for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
		int	valor1 = -1,
			valor2 = -1;
		for (int j = P - 1 + H; j < P - 1 + H + O; j++) {
			if (pai.ind[j] == filho.ind[i])
				valor1 = pai.n_vezes[j];
			if (mae.ind[j] == filho.ind[i])
				valor2 = mae.n_vezes[j];
			if (valor1 != -1 && valor2 != -1)
				break;
		}
		filho.n_vezes[i] = ceil((double)(valor1 + valor2) / 2);
	}

	default_random_engine generator;
	generator.seed(time(NULL));
	uniform_real_distribution<double> distribution(0.0, 1.0);

	if (!viavel(filho))
		corrigir(filho);

	if (distribution(generator) < prob_mutacao) {
		mutar(filho);
		if (!viavel(filho))
			corrigir(filho);
	}

	if (viavel(filho)) {
		filho.fitness = fitness(filho);
		filhos.push_back(filho);
	}



	//GERAÇÃO DO FILHO 2
	filho = mae;
	
	for (int i = 0; i < P - 1 + H + O; i++) {
		filho.ind[i] = -1;
		filho.n_vezes[i] = -1;
	}


	//CROSSOVER INDICES DE PADRÕES DE EMPACOTAMENTO
	for (int i = 0; i < floor((double)(P - 1) / 2); i++)
		filho.ind[i] = mae.ind[i];

	nao_escolhidos.clear();
	nao_escolhidos.resize(0);


	for (int i = floor((double)(P - 1) / 2); i < P - 1; i++)
		nao_escolhidos.push_back(mae.ind[i]);

	for (int i = floor((double)(P - 1) / 2); i < P - 1; i++) {
		if (find(nao_escolhidos.begin(), nao_escolhidos.end(), pai.ind[i]) != nao_escolhidos.end()) {
			filho.ind[i] = pai.ind[i];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
		else {
			default_random_engine generator;
			generator.seed(time(NULL));
			uniform_int_distribution<int> distribution(0, nao_escolhidos.size() - 1);

			int escolhido = distribution(generator);

			filho.ind[i] = nao_escolhidos[escolhido];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
	}

	//CROSSOVER INDICES DE PADRÕES DE CORTE
	for (int i = P - 1; i < P - 1 + floor((double)H / 2); i++) {
		filho.ind[i] = mae.ind[i];
	}


	for (int i = P - 1 + floor((double)H / 2); i < P - 1 + H; i++)
		nao_escolhidos.push_back(mae.ind[i]);

	for (int i = P - 1 + floor((double)H / 2); i < P - 1 + H; i++) {
		if (find(nao_escolhidos.begin(), nao_escolhidos.end(), pai.ind[i]) != nao_escolhidos.end()) {
			filho.ind[i] = pai.ind[i];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
		else {
			default_random_engine generator;
			generator.seed(time(NULL));
			uniform_int_distribution<int> distribution(0, nao_escolhidos.size() - 1);

			int escolhido = distribution(generator);

			filho.ind[i] = nao_escolhidos[escolhido];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
	}

	//CROSSOVER INDICES DE PADRÕES DE TRASPASSE
	for (int i = P - 1 + H; i <P - 1 + H + floor((double)O / 2); i++) {
		filho.ind[i] = mae.ind[i];
	}
	for (int i = P - 1 + H + floor((double)O / 2); i < P - 1 + H + O; i++)
		nao_escolhidos.push_back(mae.ind[i]);

	for (int i = P - 1 + H + floor((double)O / 2); i < P - 1 + H + O; i++) {
		if (find(nao_escolhidos.begin(), nao_escolhidos.end(), pai.ind[i]) != nao_escolhidos.end()) {
			filho.ind[i] = pai.ind[i];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
		else {
			default_random_engine generator;
			generator.seed(time(NULL));
			uniform_int_distribution<int> distribution(0, nao_escolhidos.size() - 1);

			int escolhido = distribution(generator);

			filho.ind[i] = nao_escolhidos[escolhido];
			remove(nao_escolhidos.begin(), nao_escolhidos.end(), filho.ind[i]);
			nao_escolhidos.resize(nao_escolhidos.size() - 1);
		}
	}

	//CALCULAR NOVAS QUANTIDADES PELAS MÉDIAS
	for (int i = 0; i < P - 1; i++) {
		int	valor1 = -1,
			valor2 = -1;
		for (int j = 0; j < P - 1; j++) {
			if (pai.ind[j] == filho.ind[i])
				valor1 = pai.n_vezes[j];
			if (mae.ind[j] == filho.ind[i])
				valor2 = mae.n_vezes[j];
			if (valor1 != -1 && valor2 != -1)
				break;
		}
		filho.n_vezes[i] = ceil((double)(valor1 + valor2) / 2);
	}

	for (int i = P - 1; i < P - 1 + H; i++) {
		int	valor1 = -1,
			valor2 = -1;
		for (int j = P - 1; j < P - 1 + H; j++) {
			if (pai.ind[j] == filho.ind[i])
				valor1 = pai.n_vezes[j];
			if (mae.ind[j] == filho.ind[i])
				valor2 = mae.n_vezes[j];
			if (valor1 != -1 && valor2 != -1)
				break;
		}
		filho.n_vezes[i] = ceil((double)(valor1 + valor2) / 2);
	}

	for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
		int	valor1 = -1,
			valor2 = -1;
		for (int j = P - 1 + H; j < P - 1 + H + O; j++) {
			if (pai.ind[j] == filho.ind[i])
				valor1 = pai.n_vezes[j];
			if (mae.ind[j] == filho.ind[i])
				valor2 = mae.n_vezes[j];
			if (valor1 != -1 && valor2 != -1)
				break;
		}
		filho.n_vezes[i] = ceil((double)(valor1 + valor2) / 2);
	}


	if (!viavel(filho))
		corrigir(filho);

	if (distribution(generator) < prob_mutacao) {
		mutar(filho);
		if (!viavel(filho))
			corrigir(filho);
	}


	if (viavel(filho)) {
		filho.fitness = fitness(filho);
		filhos.push_back(filho);
	}


	return filhos;
}





list<individuo> Heuristica::cruzamento(vector<individuo> Popu)
{
	default_random_engine generator;
	generator.seed(time(NULL));
	uniform_real_distribution<double> distribution(0.0, 1.0);

	list<individuo> retorno;
	for (int i = 0; i < Popu.size(); i++) {
		for (int j = 0; j < Popu.size(); j++) {
			if (i != j && distribution(generator) < prob_cruzamento) {
				list<individuo> aux = cruzar(Popu[i], Popu[j]);
				retorno.insert(retorno.end(), aux.begin(), aux.end());
			}
		}
	}

	return retorno;
}








bool Heuristica::viavel(individuo solu)
{
	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;

	//Inviabilidade por permutação
	/*vector<int> p1(P - 1), p2(H), p3(O),
				s1(P - 1), s2(H), s3(O);

	for (int i = 0; i < P - 1; i++)
		s1[i] = solu.ind[i];
	for (int i = P - 1; i < P - 1 + H; i++)
		s2[i - (P - 1)] = solu.ind[i];
	for (int i = P - 1 + H; i < P - 1 + H + O; i++)
		s3[i - (P - 1 + H)] = solu.ind[i];

	iota(p1.begin(), p1.end(), 1);
	iota(p2.begin(), p2.end(), 0);
	iota(p3.begin(), p3.end(), 0);

	if (!is_permutation(p1.begin(), p1.end(), s1.begin()))
		return false;
	if (!is_permutation(p2.begin(), p2.end(), s2.begin()))
		return false;
	if (!is_permutation(p3.begin(), p3.end(), s3.begin()))
		return false;*/


	//Iniciar as demandas atuais em 0
	for (auto &elemento : DemandasAuxiliares)
		for (auto &demand : elemento.demandas)
			demand = 0;

	for (int i = 0; i < P - 1; i++)
	{
		//Contar as demandas que o padrão atual preenche
		for (int tam = 0; tam < PackPatterns[solu.ind[i]].n_comprimentos; tam++) {
			DemandasAuxiliares[PackPatterns[solu.ind[i]].tipo].demandas[tam]
				+= solu.n_vezes[i] * PackPatterns[solu.ind[i]].tamanhos[tam];
		}
	}


	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++)
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
			if (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k])
				return false;
	
	for (int i = 0; i < P - 1; i++)
	{
		if (solu.n_vezes[i] == 0)
			continue;

		for (int gamma = 0; gamma < Gamma; gamma++)
			if (PackPatterns[solu.ind[i]].maximal(L[gamma], TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos))
				FormasQueDevemSerGeradas[gamma] += TipoVigas[PackPatterns[solu.ind[i]].tipo].n_barras * solu.n_vezes[i];
	}

	
	//Viabilidade por estoque

	vector<int> EstoqueUsado(W + V, 0);



	for (int i = P - 1; i < P - 1 + H + O; i++) {
		if (i < P - 1 + H)
			EstoqueUsado[CutPatterns[solu.ind[i]].index_barra] += solu.n_vezes[i];
		else {
			for (int v = 0; v < V; v++)
				EstoqueUsado[W + v] += solu.n_vezes[i] * SplPatterns[solu.ind[i]].tamanhos[v];
		}
	}

	for (int w = 0; w < W + V; w++) {
		if (EstoqueUsado[w] > e[w]) {
			return false;
		}
	}


	//contar as barras com tamanhos de forma geradas
	vector<int> BarrasGeradasPelaSolu(Gamma, 0);
	for (int gamma = 0; gamma < Gamma; gamma++)
		BarrasGeradasPelaSolu[gamma] = 0;


	for (int i = P - 1; i < P - 1 + H; i++) {
		for (int gamma = 0; gamma < Gamma; gamma++) {
			if(CutPatterns[solu.ind[i]].tamanhos[gamma] != 0)
				BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
		}
	}

	for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
		for (int gamma = 0; gamma < Gamma; gamma++) {
			if (SplPatterns[solu.ind[i]].barra_gerada == gamma)
				BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
		}
	}


	//Se a quantidade da solu for diferente da que deveria ser é inviável
	for (int gamma = 0; gamma < Gamma; gamma++)
		if (BarrasGeradasPelaSolu[gamma] != FormasQueDevemSerGeradas[gamma])
			return false;

	return true;
}







void Heuristica::ImprimirSolucaoEstiloCPLEX(individuo solu) {

	return;
}




void Heuristica::corrigir(individuo &solu)
{
	
	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;
	//Iniciar as demandas atuais em 0
	for (auto &elemento : DemandasAuxiliares)
		for (auto &demand : elemento.demandas)
			demand = 0;


	bool ja_encheu = false;

	for (int i = 0; i < P - 1; i++)
	{
		if (ja_encheu && solu.n_vezes[i] > 0)
			solu.n_vezes[i] = 0;
		else {

			for (int nvezes = 0; nvezes < solu.n_vezes[i]; nvezes++) {

				//Contar as demandas que o padrão atual preenche
				for (int tam = 0; tam < PackPatterns[solu.ind[i]].n_comprimentos; tam++) {
					DemandasAuxiliares[PackPatterns[solu.ind[i]].tipo].demandas[tam]
						+= PackPatterns[solu.ind[i]].tamanhos[tam];
				}
				bool entrou = false;
				for (int c = 0; c < C; c++)
					for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
						if (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k])
							entrou = true;
				if (!entrou) {
					ja_encheu = true;
					solu.n_vezes[i] = nvezes + 1;
					break;
				}
			}
		}
	}


	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++) {
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++) {
			//Se a demanda do tipo c tamanho k não é atendida
			if (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k]) {
				//Pegar o primeiro padrão que cobre ela
				int primeiro_padrao = - 1;
				for (int i = 0; i < P - 1; i++) {
					if (PackPatterns[solu.ind[i]].tipo == c && PackPatterns[solu.ind[i]].tamanhos[k] > 0) {
						primeiro_padrao = i;
						break;
					}
				}
				int n_add = ceil((double) (TipoVigas[c].demandas[k] - DemandasAuxiliares[c].demandas[k])
					/PackPatterns[solu.ind[primeiro_padrao]].tamanhos[k]);
				
				solu.n_vezes[primeiro_padrao] += n_add;
				for (int k2 = 0; k2 < TipoVigas[c].n_comprimentos; k2++)
					DemandasAuxiliares[c].demandas[k2] += PackPatterns[solu.ind[primeiro_padrao]].tamanhos[k2];
				
			}
		}
	}
	//Contar qtd necessária de barras
	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	for (int i = 0; i < P - 1; i++)
	{
		int tipo_do_padrao = PackPatterns[solu.ind[i]].tipo;
		for (int gamma = 0; gamma < Gamma; gamma++)
		{
			if (PackPatterns[solu.ind[i]].maximal(L[gamma], TipoVigas[tipo_do_padrao].comprimentos))
				FormasQueDevemSerGeradas[gamma] += TipoVigas[tipo_do_padrao].n_barras * solu.n_vezes[i];
		}
	}

	vector<int> EstoqueUsado(W+V, 0);
	//Corrigir por estoque de barras

	for (int i = P - 1; i < P - 1 + H + O; i++) {
		if (i < P - 1 + H)
			EstoqueUsado[CutPatterns[solu.ind[i]].index_barra] += solu.n_vezes[i];
		else {
			for (int v = 0; v < V; v++)
				EstoqueUsado[W + v] += solu.n_vezes[i] * SplPatterns[solu.ind[i]].tamanhos[v];
		}
	}


	for (int w = 0; w < W + V; w++) {
		//se desobedece o estoque
		if (EstoqueUsado[w] > e[w]) {
			//Calcular qtde que deve ser retirada
			int retirar = EstoqueUsado[w] - e[w];
			//Para todo padrão na ordem que está na solução
			for (int i = P - 1; i < P - 1 + H; i++) {
				//se o padrão corta a barra w e ele está na solução
				if (CutPatterns[solu.ind[i]].index_barra == w && solu.n_vezes[i] > 0) {
					/*qtde que vai ser removidade do padrão da solução é o menor valor entre o que
					tem que ser removido e o número de vezes que o padrão é usado
					*/
					int remover = min(solu.n_vezes[i], retirar);
					//atualiza o número de vezes que o padrão tá na solução
					solu.n_vezes[i] -= remover;
					//atualiza o estoque usado da barra w
					EstoqueUsado[w] -= remover;
					//atualiza a qtde que ainda falta ser retirada
					retirar -= remover;
					//se já removeu o que tinha q remover para de percorrer os padrões
					if (retirar <= 0)
						break;
				}
			}
			//Se mesmo depois de percorrer os padrões ainda houver barras a serem removidas, percorra os padrões de traspasse
			if (retirar > 0) {
				for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
					//Se o padrão de traspasse está na solução
					if (solu.n_vezes[i] > 0) {
						/*qtde a ser removida é o minimo entre o número de vezes que ele está na solução e o piso
						do número que preciso retirar sobre a qtde de barras w que o padrão usa
						*/
						int remover = min(solu.n_vezes[i],
							int(floor((double) retirar / SplPatterns[solu.ind[i]].tamanhos[w - W])));
						solu.n_vezes[i] -= remover;

						for (int v = 0; v < V; v++)
							EstoqueUsado[w] -= SplPatterns[solu.ind[i]].tamanhos[v - W] * remover;

						retirar -= remover;
						if (retirar <= 0)
							break;
					}
				}
			}
		}
	}
	//Tudo OK até aqui
	
	vector<int> BarrasGeradas(Gamma, 0);
	//Contar barras que estão sendo geradas
	for (int gamma = 0; gamma < Gamma; gamma++) {
		for (int i = P - 1; i < P - 1 + H + O; i++) {
			if (i < P - 1 + H)
				BarrasGeradas[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
			else {
				if(SplPatterns[solu.ind[i]].barra_gerada == gamma)
					BarrasGeradas[gamma] += solu.n_vezes[i];
			}
		}
	}
	
	for (int gamma = 0; gamma < Gamma; gamma++) {

		if (BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {

			//Se tem em excesso, remover dos padrões que só tem o tipo gamma
			for (int i = P - 1; i < P - 1 + H; i++) {
				if (solu.n_vezes[i] > 0) {
					//só gera formas do tipo gamma?
					bool gera_gamma = true;
					for (int tam = 0; tam < Gamma; tam++) {
						if (tam != gamma && CutPatterns[solu.ind[i]].tamanhos[tam] > 0) {
							gera_gamma = false;
							break;
						}
					}
					if (gera_gamma && CutPatterns[solu.ind[i]].tamanhos[gamma] == 0)
						gera_gamma = false;


					if (gera_gamma) {

						int remover = min(solu.n_vezes[i], (int)ceil((double)(BarrasGeradas[gamma] - FormasQueDevemSerGeradas[gamma])
							/ CutPatterns[solu.ind[i]].tamanhos[gamma]));

						solu.n_vezes[i] -= remover;
						BarrasGeradas[gamma] -= remover * CutPatterns[solu.ind[i]].tamanhos[gamma];
						EstoqueUsado[CutPatterns[solu.ind[i]].index_barra] -= remover;

						
					}
					if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
						break;
				}

			}
		

			//Se continuar, remova dos traspasses
			if (BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {
				for (int i = P - 1 + H; i < P - 1 + H + O; i++) {

					if (SplPatterns[solu.ind[i]].barra_gerada == gamma && solu.n_vezes[i] > 0) {
						int remover = min(solu.n_vezes[i], BarrasGeradas[gamma] - FormasQueDevemSerGeradas[gamma]);
						solu.n_vezes[i] -= remover;
						BarrasGeradas[gamma] -= remover;

						for (int v = 0; v < V; v++)
							EstoqueUsado[W + v] -= remover * SplPatterns[solu.ind[i]].tamanhos[v];	
					}
					if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
						break;
				}
			}
		}



		if (BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma]) {

			//Se ta em falta, adiciona
			for (int i = P - 1; i < P - 1 + H; i++) {
				//Se gera só gamma
				bool gera_gamma = true;
				for (int tam = 0; tam < Gamma; tam++) {
					if (tam != gamma && CutPatterns[solu.ind[i]].tamanhos[tam] > 0) {
						gera_gamma = false;
						break;
					}
				}
				if (gera_gamma && CutPatterns[solu.ind[i]].tamanhos[gamma] == 0)
					gera_gamma = false;


				//se o padrão cobre uma barra de tamanho gamma
				if (gera_gamma) {
					//enquanto ainda é menor e o padrão não interfere no estoque da barra
					int adicionar = min((int)(floor((double)(FormasQueDevemSerGeradas[gamma] - BarrasGeradas[gamma]) / CutPatterns[solu.ind[i]].tamanhos[gamma])),
						e[CutPatterns[solu.ind[i]].index_barra] - EstoqueUsado[CutPatterns[solu.ind[i]].index_barra]);

					EstoqueUsado[CutPatterns[solu.ind[i]].index_barra] += adicionar;
					solu.n_vezes[i] += adicionar;

					BarrasGeradas[gamma] += adicionar * CutPatterns[solu.ind[i]].tamanhos[gamma];

					if (BarrasGeradas[gamma] == FormasQueDevemSerGeradas[gamma])
						break;
				}
			}

			//Se ainda não supriu adiciona splicing feito doido
			if (BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma]) {
				for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
					if (SplPatterns[solu.ind[i]].barra_gerada == gamma) {
						bool pode_ser_adicionado = true;
						while (pode_ser_adicionado  && BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma]) {

							for (int v = 0; v < V; v++) {
								if (EstoqueUsado[W + v] + SplPatterns[solu.ind[i]].tamanhos[v] > e[W + v]) {
									pode_ser_adicionado = false;
									break;
								}
							}
							if (pode_ser_adicionado) {
								solu.n_vezes[i]++;
								BarrasGeradas[gamma]++;
								for (int v = 0; v < V; v++)
									EstoqueUsado[W + v] += SplPatterns[solu.ind[i]].tamanhos[v];
							}
						}
						if (BarrasGeradas[gamma] == FormasQueDevemSerGeradas[gamma])
							break;
					}
				}
			}

		}
	}

	if (!viavel(solu)) {
		cout << "Filho que nao conseguiu corrigir detectado!" << endl;
	}

}


int Heuristica::qtde_adicionavel(Padrao_Traspasse Padrao, vector<int> EstoqueUsado, list<int> usados) {

	vector<int> qtdes(usados.size(), 0);

	int i = 0;
	for (auto w : usados) {
		qtdes[i] = floor((e[w] - EstoqueUsado[w])/Padrao.tamanhos[w - W]);
		i++;
	}

	return *min_element(qtdes.begin(), qtdes.end());
}





void Heuristica::ImprimirVetorSolu(individuo solu)
{

	for (int i = 0; i < P - 1; i++) {
		if(solu.n_vezes[i] != 0)
		cout << PackPatterns[solu.ind[i]].id << "," << solu.n_vezes[i] << " ";
	}
	cout << "//" << endl;
	for (int i = P - 1; i < P - 1 + H; i++) {
		if (solu.n_vezes[i] != 0)
			cout << CutPatterns[solu.ind[i]].index_pat << "," << solu.n_vezes[i] << " ";
	}
	cout << "//" << endl;
	for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
		if (solu.n_vezes[i] != 0)
			cout << SplPatterns[solu.ind[i]].id << "," << solu.n_vezes[i] << " ";
	}

	return;





	for (int i = 0; i < P - 1 + H + O; i++) {
		if (i < P - 1)
			cout << PackPatterns[solu.ind[i]].id << "," << solu.n_vezes[i] << " ";
		else {
			if (i < P - 1 + H)
				cout << CutPatterns[solu.ind[i]].index_pat << "," << solu.n_vezes[i] << " ";
			else
				cout << SplPatterns[solu.ind[i]].id << "," << solu.n_vezes[i] << " ";
		}
	}
}





individuo Heuristica::GerarSoluGRASP() {
	vector<int>	solu1(P - 1 + H + O),	//índices dos padrões
				solu2(P - 1 + H + O, 0);	//quantidade dos padrões

	//Inicializando os índices

	iota(solu1.begin(), solu1.begin() + P - 1, 1);
	iota(solu1.begin() + P - 1, solu1.begin() + P - 1 + H, 0);
	iota(solu1.begin() + P - 1 + H, solu1.end(), 0);
	

	//Fazendo permutações aleatórias dos índices
	random_shuffle(solu1.begin(), solu1.begin() + P - 1);
	random_shuffle(solu1.begin() + P - 1, solu1.begin() + P - 1 + H);
	random_shuffle(solu1.begin() + P - 1 + H, solu1.end());
	

	//prenchendo os padrões de empacotamento até a demanda ser atendida
	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;

	//Iniciar as demandas atuais em 0
	for (auto &elemento : DemandasAuxiliares)
		for (auto &demand : elemento.demandas)
			demand = 0;

	for (int i = 0; i < P - 1; i++) {
		int tipo_atual = PackPatterns[solu1[i]].tipo;
		
		bool demanda_tipo_atendida = true;
		list<int> SemDemanda_PadraoCobre;
		vector<int> qtde_necessaria(TipoVigas[tipo_atual].n_comprimentos, 0);
		for (int k = 0; k < TipoVigas[tipo_atual].n_comprimentos; k++) {
			if (DemandasAuxiliares[tipo_atual].demandas[k] < TipoVigas[tipo_atual].demandas[k]) {
				demanda_tipo_atendida = false;
				if(PackPatterns[solu1[i]].tamanhos[k] > 0)
					SemDemanda_PadraoCobre.push_back(k);
				qtde_necessaria[k] = TipoVigas[tipo_atual].demandas[k] - DemandasAuxiliares[tipo_atual].demandas[k];
			}
		}
		if (demanda_tipo_atendida)
			continue;


		if (!demanda_tipo_atendida) {
			int maior_necessaria = -1;
			//teste se o padrão cobre o tamanho necessario
			//pega o tamanho mais necessário que o padrão cobre

			for (auto k_aux : SemDemanda_PadraoCobre) {
				int n_vezes_uso_padrao_atual = ceil((double)(qtde_necessaria[k_aux] / PackPatterns[solu1[i]].tamanhos[k_aux]));

				solu2[i] += n_vezes_uso_padrao_atual;

				for (int k = 0; k < TipoVigas[tipo_atual].n_comprimentos; k++) {
					DemandasAuxiliares[tipo_atual].demandas[k] +=
						n_vezes_uso_padrao_atual*PackPatterns[solu1[i]].tamanhos[k];
					qtde_necessaria[k] = TipoVigas[tipo_atual].demandas[k] - DemandasAuxiliares[tipo_atual].demandas[k];
				}

			}

		}
	}

	//Populando a quantidade de barras que devem ser produzidas
	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	for (int i = 0; i < P - 1; i++) {
		for (int gamma = 0; gamma < Gamma; gamma++)
			if (PackPatterns[solu1[i]].maximal(L[gamma], TipoVigas[PackPatterns[solu1[i]].tipo].comprimentos))
				FormasQueDevemSerGeradas[gamma] += TipoVigas[PackPatterns[solu1[i]].tipo].n_barras*solu2[i];
	}





	vector<int> FormasGeradas(Gamma, 0);
	vector<int> EstoqueDisponivel = e;

	//cout << "\n\n\n";
	for (int gamma = 0; gamma < Gamma; gamma++) {
		int necessario = FormasQueDevemSerGeradas[gamma];

		for (int i = P - 1; i < P - 1 + H; i++) {
			if (CutPatterns[solu1[i]].tamanhos[gamma] > 0 && necessario > 0) {

				int adicionado = floor(min(EstoqueDisponivel[CutPatterns[solu1[i]].index_barra], necessario) / CutPatterns[solu1[i]].tamanhos[gamma]);
				FormasGeradas[gamma] += CutPatterns[solu1[i]].tamanhos[gamma] * adicionado;
				solu2[i] += adicionado;
				EstoqueDisponivel[CutPatterns[solu1[i]].index_barra] -= adicionado;
				necessario -= CutPatterns[solu1[i]].tamanhos[gamma] * adicionado;
			}
		}
		//se já conseguiu preencher com padrãoes de corte não precisa tentar com os traspasses
		if (necessario == 0)
			continue;


		for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
			if (SplPatterns[solu1[i]].barra_gerada == gamma) {
				bool pode_ser_adicionado = true;
				while (pode_ser_adicionado  && necessario > 0) {

					for (int v = 0; v < V; v++) {
						if (EstoqueDisponivel[W + v] - SplPatterns[solu1[i]].tamanhos[v] < 0) {
							pode_ser_adicionado = false;
							break;
						}
					}
					if (pode_ser_adicionado) {
						solu2[i]++;
						FormasGeradas[gamma]++;
						necessario--;
						for (int v = 0; v < V; v++) {
							EstoqueDisponivel[W + v] -= SplPatterns[solu1[i]].tamanhos[v];
						}
					}
				}
				if (FormasGeradas[gamma] == FormasQueDevemSerGeradas[gamma])
					break;
			}
		}



		
	}

	//Se com os padrões de corte não supriu o necessário de barras adicione splicing até suprir

	vector<int> EstoqueDisponivel_Aux;

	for (int gamma = 0; gamma < Gamma; gamma++) {
		int necessario = FormasQueDevemSerGeradas[gamma] - FormasGeradas[gamma];
		for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
			if (SplPatterns[solu1[i]].barra_gerada == gamma && necessario > 0) {
				EstoqueDisponivel_Aux = EstoqueDisponivel;

				int adicionado = necessario;

				bool pular_padrao = false;
				for (int v = 0; v < V; v++) {
					EstoqueDisponivel_Aux[W + v] -= SplPatterns[solu1[i]].tamanhos[v] * adicionado;
					if (EstoqueDisponivel_Aux[W + v] < 0) {
						if(adicionado > floor((double)EstoqueDisponivel[W + v] / SplPatterns[solu1[i]].tamanhos[v]))
							adicionado = floor((double)EstoqueDisponivel[W + v] / SplPatterns[solu1[i]].tamanhos[v]);
					}
				}
				
				for (int v = 0; v < V; v++) 
					EstoqueDisponivel[W + v] -= SplPatterns[solu1[i]].tamanhos[v] * adicionado;


				solu2[i] += adicionado;
				necessario -= adicionado;
				FormasGeradas[gamma] += adicionado;
			}
		}
	}

	individuo solucao;
	solucao.ind = solu1;
	solucao.n_vezes = solu2;
	return solucao;
}





void Heuristica::selecao(vector<individuo> &Popu)	//Selecao por ranking
{

	sort(Popu.begin(), Popu.end(), [](individuo i1, individuo i2) {return i1.fitness < i2.fitness; });
	vector<individuo> Auxiliar = Popu;
	Popu.clear();
	Popu.resize(0);

	for(int j = 0; j < 1; j++)
		list_tabu.push_back(Auxiliar[j]);

	if (list_tabu.size() > 20) {
		for (int j = 0; j < 20 - list_tabu.size(); j++)
			list_tabu.pop_back();
	}
	
	
	//Elistimo
	Popu.push_back(Auxiliar[0]);

	//ranking
	int i = 1,
		contador = 1;
	while (contador < TamanhoDaPopulacao) {
		bool tabu = false;

		for (auto tabuzinho : list_tabu) {
			if (Auxiliar[i].comparacao(tabuzinho) > 0.75) {
				tabu = true;
				break;
			}
		}
		if (!tabu) {
			Popu.push_back(Auxiliar[i]);
			contador++;
		}

		i++;
	}
}






void Heuristica::mutar(individuo & solu){
	//return;
	default_random_engine generator;
	generator.seed(time(NULL));
	uniform_int_distribution<int> bloco1(0, P - 2);
	//indices que serao trocados
	int id_pack1, id_pack2;

	vector<int> diferentes_de_zero;
	for (int i = 0; i < P - 1; i++)
		if (solu.n_vezes[i] > 0)
			diferentes_de_zero.push_back(i);

	uniform_int_distribution<int> bloco1_n0(0, diferentes_de_zero.size() - 1);
	int aux_id = bloco1_n0(generator);
	id_pack1 = diferentes_de_zero[aux_id];
	do {
		id_pack2 = bloco1(generator);
	} while (id_pack2 == id_pack1 || PackPatterns[solu.ind[id_pack1]].tipo != PackPatterns[solu.ind[id_pack2]].tipo);

	int aux = solu.ind[id_pack1];
	solu.ind[id_pack1] = solu.ind[id_pack2];
	solu.ind[id_pack2] = aux;




	int id_cut1, id_cut2;
	diferentes_de_zero.clear();
	for (int i = P - 1; i < P - 1 + H; i++)
		if (solu.n_vezes[i] > 0)
			diferentes_de_zero.push_back(i);

	uniform_int_distribution<int> bloco2_n0(0, diferentes_de_zero.size() - 1);
	uniform_int_distribution<int> bloco2(P - 1, P - 2 + H);
	aux_id = bloco2_n0(generator);
	id_cut1 = diferentes_de_zero[aux_id];
	solu.n_vezes[id_cut1] = 0;

	do {
		id_cut2 = bloco2(generator);
	} while (id_cut1 == id_cut2);

	aux = solu.ind[id_cut1];
	solu.ind[id_cut1] = solu.ind[id_cut2];
	solu.ind[id_cut2] = aux;
}





void Heuristica::restart_populacao(vector<individuo> &Populacao) {
	individuo Elite = Populacao[0];
	Populacao.clear();
	Populacao.resize(0);
	Populacao.push_back(Elite);
	list_tabu.clear();

	for (int i = 0; i < P + H + O; i++)
	{
		individuo solucao = GerarSoluGRASP();

		//corrigir(solucao);

		if (!viavel(solucao)) {
			cout << "Erro! Solucao gulosa gerada inviavel!" << endl;
			//system("pause");
			//exit(0);
		}
		else {
			solucao.fitness = fitness(solucao);
			Populacao.push_back(solucao);
		}
	}
}






void Heuristica::funcaoteste() {
	srand(time(NULL));
	definir_parametros();
	vector<individuo> Populacao;

	cout << "Gerando populacao inicial\n";

	auto start = chrono::system_clock::now();
	
	for (int i = 0; i < P + H + O; i++)
	{
		individuo solucao = GerarSoluGRASP();

		//corrigir(solucao);
		
		if (!viavel(solucao)) {
			cout << "Erro! Solucao gulosa gerada inviavel!" << endl;
			//system("pause");
			//exit(0);
		}
		else {
			solucao.fitness = fitness(solucao);
			Populacao.push_back(solucao);
		}	
	}
	

	selecao(Populacao);
	ImprimirVetorSolu(Populacao[0]);
	cout << "Populacao inicial gerada e elitizada\n";

	for (auto solucao : Populacao)
		cout << "FO = " << solucao.fitness << endl;


	double	fit_atual,
			fit_antigo;
	int contador_fit = 0;


	for(int i = 0; i < NGeracoes; ++i){
		fit_antigo = Populacao[0].fitness;

		list<individuo> Offspring = cruzamento(Populacao);
		Populacao.insert(Populacao.end(), Offspring.begin(), Offspring.end());
		selecao(Populacao);

		fit_atual = Populacao[0].fitness;



		cout << "\n\n Geracao" << i << endl;
		for (auto solucao : Populacao)
			cout << "FO = " << solucao.fitness << endl;
		ImprimirVetorSolu(Populacao[0]);

		if(viavel(Populacao[0]))
			cout << "Viavel!!";

		
		if(i%2 != 0)
			prob_mutacao += taxa_aumento_mut*prob_mutacao;

		if (fit_atual == fit_antigo) {
			contador_fit++;
			if (contador_fit > 10) {
				contador_fit = 0;
				//restart_populacao(Populacao);
				//prob_mutacao = 0.05;
				break;
			}
		}
		else {
			contador_fit = 0;
		}
	}
	

	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "\t Tempo Total " << elapsed_seconds.count() << "s" << endl;

	ImprimirVetorSolu(Populacao[0]);
}