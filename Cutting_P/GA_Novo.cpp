#include "GA_Novo.h"

//CALCULA FUNÇÃO OBJETIVO
double GA_Novo::fitness(individuo solu){

	int tamanho_gene = solu.ind.size();
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


	for (int i = 0; i < tamanho_gene; i++)
	{
		if (solu.ind[i] <= P - 1) {
			int tipo_da_forma = Gamma_Associado[solu.ind[i]];
			
			for (int vezes = 0; vezes < solu.n_vezes[i]; vezes++) {
				std::sort(UtilFormas.begin(), UtilFormas.end());
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
		else if (solu.ind[i] <= P - 1 + H) {
			int indice_pat = solu.ind[i] - P;
			if (CutPatterns[indice_pat].index_barra < W && !CutPatterns[indice_pat].Gera_LeftOvers(Gamma, V))
				Sobra1 += solu.n_vezes[i] * (b[CutPatterns[indice_pat].index_barra] - CutPatterns[indice_pat].cap);

			if (CutPatterns[indice_pat].index_barra < W && CutPatterns[indice_pat].Gera_LeftOvers(Gamma, V))
				Sobra2 += solu.n_vezes[i] * (b[CutPatterns[indice_pat].index_barra] - CutPatterns[indice_pat].cap);

			if (CutPatterns[indice_pat].index_barra >= W)
				Sobra3 += solu.n_vezes[i] * (b[CutPatterns[indice_pat].index_barra] - CutPatterns[indice_pat].cap);
		}
		else if(solu.ind[i] <= P - 1 + H + O) {
			int indice_pat = solu.ind[i] - P - H;
			Sobra3 += solu.n_vezes[i] * (SplPatterns[indice_pat].folga);
		}
		
	}
	Utilizacao valor = *max_element(UtilFormas.begin(), UtilFormas.end());
	Makespan = valor.util;

	return Makespan + Sobra1 + Parameter_alpha_1*Sobra2 + Parameter_alpha_2*Sobra3;
}

//GERA SOLUÇÃO ALEATÓRIA
GA_Novo::individuo GA_Novo::GerarSoluGRASP() {
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
				if (PackPatterns[solu1[i]].tamanhos[k] > 0)
					SemDemanda_PadraoCobre.push_back(k);
				qtde_necessaria[k] = TipoVigas[tipo_atual].demandas[k] - DemandasAuxiliares[tipo_atual].demandas[k];
			}
		}
		if (demanda_tipo_atendida || SemDemanda_PadraoCobre.size() == 0)
			continue;


		if (!demanda_tipo_atendida) {
			int maior_necessaria = -1;
			//teste se o padrão cobre o tamanho necessario
			//pega o tamanho mais necessário que o padrão cobre

			for (auto k_aux : SemDemanda_PadraoCobre) {
				int n_vezes_uso_padrao_atual = ceil((double)qtde_necessaria[k_aux] / PackPatterns[solu1[i]].tamanhos[k_aux]);

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

	//std::cout << "\n\n\n";
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
						if (adicionado > floor((double)EstoqueDisponivel[W + v] / SplPatterns[solu1[i]].tamanhos[v]))
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
	solucao.ind = vector<int>(0, 0);
	solucao.n_vezes = vector<int>(0, 0);

	for (int i = 0; i < P - 1 + H + O; i++) {
		if (solu2[i] > 0) {
			solucao.n_vezes.push_back(solu2[i]);
			if (i < P - 1)
				solucao.ind.push_back(solu1[i]);	
			else if (i < P - 1 + H)
				solucao.ind.push_back(solu1[i] + P);
			else if(i < P - 1 + H + O)
				solucao.ind.push_back(solu1[i] + P + H);
		}
	}
	return solucao;
}

//TESTA VIABILIDADE DE UM INDIVÍDUO
bool GA_Novo::viavel(individuo solu)
{
	int n_genes = solu.n_vezes.size();

	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;

	//Iniciar as demandas atuais em 0
	for (auto &elemento : DemandasAuxiliares)
		for (auto &demand : elemento.demandas)
			demand = 0;

	for (int i = 0; i < n_genes; i++)
	{
		//Se caracteriza um padrao de empacotamento
		if (solu.ind[i] <= P - 1) {
			//Contar as demandas que o padrão atual preenche
			for (int tam = 0; tam < PackPatterns[solu.ind[i]].n_comprimentos; tam++) {
				DemandasAuxiliares[PackPatterns[solu.ind[i]].tipo].demandas[tam]
					+= solu.n_vezes[i] * PackPatterns[solu.ind[i]].tamanhos[tam];
			}
		}
	}


	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++)
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
			if (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k])
				return false;

	//Contar barras que devem ser geradas
	for (int i = 0; i < n_genes; i++)
	{
		//Se caracteriza um padrao de empacotamento
		if (solu.ind[i] <= P - 1) {
			FormasQueDevemSerGeradas[Gamma_Associado[solu.ind[i]]] += 
				TipoVigas[PackPatterns[solu.ind[i]].tipo].n_barras * solu.n_vezes[i];
		}
	}


	//Viabilidade por estoque
	vector<int> EstoqueUsado(W + V, 0);
		
	for (int i = 0; i < n_genes; i++) {
		//se representa um padrao de corte
		if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
			int indice_pat = solu.ind[i] - P;
			EstoqueUsado[CutPatterns[indice_pat].index_barra] += solu.n_vezes[i];
		}
		//se representa um padrao de splicing
		else if(solu.ind[i] >= P + H && solu.ind[i] <= P - 1 + H + O){
			int indice_pat = solu.ind[i] - P - H;
			for (int v = 0; v < V; v++)
				EstoqueUsado[W + v] += solu.n_vezes[i] * SplPatterns[indice_pat].tamanhos[v];
		}
	}

	for (int w = 0; w < W + V; w++) {
		if (EstoqueUsado[w] > e[w]) {
			return false;
		}
	}


	//contar as barras com tamanhos de forma geradas
	vector<int> BarrasGeradasPelaSolu(Gamma, 0);

	for (int i = 0; i < n_genes; i++) {
		//se representa um padrao de corte
		if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
			int indice_pat = solu.ind[i] - P;
			for (int gamma = 0; gamma < Gamma; gamma++) {
				if (CutPatterns[indice_pat].tamanhos[gamma] != 0)
					BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[indice_pat].tamanhos[gamma];
			}
		}
		//se representa um padrao de splicing
		else if (solu.ind[i] >= P + H && solu.ind[i] <= P - 1 + H + O) {
			int indice_pat = solu.ind[i] - P - H;
			for (int gamma = 0; gamma < Gamma; gamma++) {
				if (SplPatterns[indice_pat].barra_gerada == gamma)
					BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i];
			}
		}
	}

	//Se a quantidade da solu for diferente da que deveria ser é inviável
	for (int gamma = 0; gamma < Gamma; gamma++)
		if (BarrasGeradasPelaSolu[gamma] != FormasQueDevemSerGeradas[gamma])
			return false;

	return true;
}

//ALGORITMO GENÉTICO
void GA_Novo::funcao_teste()
{
	
	definir_parametros();
	vector<individuo> Populacao(0);

	auto start = chrono::system_clock::now();

	for (int i = 0; i < 100*(P + O + H); i++)
	{
		individuo solucao = GerarSoluGRASP();

		if (!viavel(solucao))
			corrigir(solucao);

		if (viavel(solucao)) {
			solucao.fitness = fitness(solucao);
			Populacao.push_back(solucao);
		}
	}
	selecao(Populacao);

	//mutar(Populacao[0]);

	double fit_antiga = Populacao[0].fitness;
	int sem_melhora = 0,
		n_restarts = 0;

	for (int i = 0; i < NGeracoes; i++) {
		individuo filho = cruzamento_unico(Populacao);
		if (viavel(filho)) {
			filho.fitness = fitness(filho);
			Populacao.push_back(filho);
		}
		selecao(Populacao);

		if (fit_antiga == Populacao[0].fitness)
			sem_melhora++;
		else
			sem_melhora = 0;

		fit_antiga = Populacao[0].fitness;


		/*if (i % 100 == 0) {
			cout << "Geracao" << i << ": " << Populacao[0].fitness << endl;
		}*/

		if (sem_melhora > 10000) {
			//cout << "Restart" << endl;
			Restart(Populacao);
			sem_melhora = 0;
			n_restarts++;
		}
		if (n_restarts == 5)
			break;
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "\t Tempo Total " << elapsed_seconds.count() << "s" << endl;
	ImprimirArquivo(Populacao[0], elapsed_seconds.count());

}

//IMPRIME REPRESENTAÇÃO DE SOLUÇÃO
void GA_Novo::ImprimirVetorSolu(individuo solu)
{
	cout << endl;
	cout << endl;
	for(int i = 0; i < solu.ind.size(); i ++)
	{
		cout << solu.ind[i] << "," << solu.n_vezes[i] << "  ";
	}
	cout << endl;
	cout << endl;
}

//SELECIONA ELEMENTOS ÚNICOS BASEADOS NA FUNÇÃO OBJETIVO
void GA_Novo::selecao(vector<individuo> &Popu)	//Selecao por ranking
{

	//Orderna populacao por fitness
	std::sort(Popu.begin(), Popu.end(), [](individuo i1, individuo i2) {return i1.fitness < i2.fitness; });
	//Salva populacao antiga
	vector<individuo> Auxiliar = Popu;
	//Limpar populacao para adicionar os individuos selecionados
	Popu.clear();

	int i = 0;
	int contador = 0;
	//adiciona sem repetição
	while (i < TamanhoDaPopulacao) {
		bool ja_entrou = false;

		for (auto elemento : Popu)
			if (elemento.igual(Auxiliar[contador]))
				ja_entrou = true;

		if (!ja_entrou) {
			Popu.push_back(Auxiliar[contador]);
			i++;
		}
		contador++;
	}
}

//GERA UM FILHO ENTRE DOIS PAIS
GA_Novo::individuo GA_Novo::cruzar(individuo pai, individuo mae)
{
	int n_pai = pai.ind.size(),
		n_mae = mae.ind.size();

	individuo filho;

	uniform_real_distribution<double> distribuicao(0.0, 1.0);

	for (int i = 0; i < n_pai; i++) {
		int ind_pai = pai.ind[i];
		int valor = pai.n_vezes[i];


		for (int j = 0; j < n_mae; j++) {
			if (mae.ind[j] == ind_pai) {
				valor = ceil((double)(valor + mae.n_vezes[j]) / 2);
				break;
			}
		}
		filho.ind.push_back(ind_pai);
		if(distribuicao(generator) > prob_mutacao)
			filho.n_vezes.push_back(valor);
		else
			filho.n_vezes.push_back(0);
	}

	for (int i = 0; i < n_mae; i++) {
		int ind_mae = mae.ind[i];
		int valor = mae.n_vezes[i];
		bool unico_na_mae = true;

		for (int j = 0; j < n_pai; j++) {
			if (pai.ind[j] == ind_mae) {
				unico_na_mae = false;
			}
		}
		if (unico_na_mae) {
			filho.ind.push_back(ind_mae);
			if (distribuicao(generator) > prob_mutacao)
				filho.n_vezes.push_back(0);
			else
				filho.n_vezes.push_back(valor);
		}
	}

	
	if (distribuicao(generator) < prob_mutacao)
		mutar(filho);

	if (!viavel(filho))
		corrigir(filho);


	return filho;
}

void GA_Novo::corrigir(individuo &solu) {
	int n_genes = solu.ind.size();

	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;
	for (auto &elemento : DemandasAuxiliares)
		for (auto &demand : elemento.demandas)
			demand = 0;

	//retirar excesso de vigas
	bool ja_encheu = false;
	for (int i = 0; i < n_genes; i++) {
		//se representa padrão de empacotamento
		if (solu.ind[i] <= P - 1) {
			//se ja atingiu a demanda
			if (ja_encheu)
				solu.n_vezes[i] = 0;
			else {
				for (int nvezes = 0; nvezes < solu.n_vezes[i]; nvezes++) {
					//Contar as demandas que o padrão atual preenche
					for (int tam = 0; tam < PackPatterns[solu.ind[i]].n_comprimentos; tam++) {
						DemandasAuxiliares[PackPatterns[solu.ind[i]].tipo].demandas[tam]
							+= PackPatterns[solu.ind[i]].tamanhos[tam];
					}
					bool n_atendeu_demanda = false;
					for (int c = 0; c < C; c++) {
						for (int k = 0; k < TipoVigas[c].n_comprimentos; k++) {
							if (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k])
								n_atendeu_demanda = true;
						}
					}
					if (!n_atendeu_demanda) {
						ja_encheu = true;
						solu.n_vezes[i] = nvezes + 1;
						break;
					}
				}
			}
		}
	}


	//corrigir por falta de vigas
	for (int c = 0; c < C; c++) {
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++) {
			//Se a demanda do tipo c tamanho k não é atendida
			while (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k]) {
				//Pegar o primeiro padrão que cobre ela
				int padrao_escolhido = -1;
				uniform_real_distribution<double> distribution(0.0, 1.0);

				for (int i = 0; i < n_genes; i++) {
					if (solu.ind[i] <= P - 1) {
						if (PackPatterns[solu.ind[i]].tipo == c && PackPatterns[solu.ind[i]].tamanhos[k] > 0 && distribution(generator) < 0.5) {
							padrao_escolhido = i;
							break;
						}
					}
				}
				if (padrao_escolhido != -1) {

					int n_add = ceil((double)(TipoVigas[c].demandas[k] - DemandasAuxiliares[c].demandas[k])
						/ PackPatterns[solu.ind[padrao_escolhido]].tamanhos[k]);

					solu.n_vezes[padrao_escolhido] += n_add;

					for (int k2 = 0; k2 < TipoVigas[c].n_comprimentos; k2++) {
						DemandasAuxiliares[c].demandas[k2] += n_add*PackPatterns[solu.ind[padrao_escolhido]].tamanhos[k2];
					}
				}
			}
		}
	}


	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	for (int i = 0; i < n_genes; i++)
	{
		//se representa um padrão de empacotamento
		if (solu.ind[i] <= P - 1) {
			int tipo_do_padrao = PackPatterns[solu.ind[i]].tipo;
			FormasQueDevemSerGeradas[Gamma_Associado[solu.ind[i]]] += TipoVigas[tipo_do_padrao].n_barras * solu.n_vezes[i];
		}
	}


	vector<int> EstoqueUsado(W + V, 0);
	//Corrigir por estoque de barras

	//Calcular estoque usado
	for (int i = 0; i < n_genes; i++) {
		if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
			int indice_real = solu.ind[i] - P;
			EstoqueUsado[CutPatterns[indice_real].index_barra] += solu.n_vezes[i];
		}
		else if (solu.ind[i] >= P) {
			int indice_real = solu.ind[i] - P - H;
			for (int v = 0; v < V; v++) {
				EstoqueUsado[W + v] += solu.n_vezes[i] * SplPatterns[indice_real].tamanhos[v];
			}
		}
	}
	
	for (int w = 0; w < W + V; w++) {
		//se desobedece o estoque
		if (EstoqueUsado[w] > e[w]) {
			//Calcular qtde que deve ser retirada
			int retirar = EstoqueUsado[w] - e[w];
			//Para todo padrão na ordem que está na solução
			for (int i = 0; i < n_genes; i++) {
				if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
					int ind_real = solu.ind[i] - P;
					//se o padrão corta a barra w e ele está na solução
					if (CutPatterns[ind_real].index_barra == w && solu.n_vezes[i] > 0) {
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
			}
			//Se mesmo depois de percorrer os padrões ainda houver barras a serem removidas, percorra os padrões de traspasse
			if (retirar > 0) {
				for (int i = 0; i < n_genes; i++) {
					if (solu.ind[i] >= P + H && solu.ind[i] <= P - 1 + H + O) {
						int ind_real = solu.ind[i] - P - H;
						//Se o padrão de traspasse está na solução
						if (solu.n_vezes[i] > 0) {
							/*qtde a ser removida é o minimo entre o número de vezes que ele está na solução e o piso
							do número que preciso retirar sobre a qtde de barras w que o padrão usa
							*/
							int remover = min(solu.n_vezes[i],
								int(floor((double)retirar / SplPatterns[ind_real].tamanhos[w - W])));
							solu.n_vezes[i] -= remover;

							for (int v = 0; v < V; v++) {
								EstoqueUsado[w] -= SplPatterns[ind_real].tamanhos[v - W] * remover;
							}
							retirar -= remover;
							if (retirar <= 0)
								break;
						}
					}
				}
			}
		}
	}

	//Contar o número de barras que está sendo gerado
	vector<int> BarrasGeradas(Gamma, 0);
	for (int gamma = 0; gamma < Gamma; gamma++) {
		for (int i = 0; i < n_genes; i++) {
			if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
				int ind_real = solu.ind[i] - P;
				BarrasGeradas[gamma] += solu.n_vezes[i] * CutPatterns[ind_real].tamanhos[gamma];
			}
			else if(solu.ind[i] >= P){
				int ind_real = solu.ind[i] - P - H;
				if (SplPatterns[ind_real].barra_gerada == gamma)
					BarrasGeradas[gamma] += solu.n_vezes[i];
			}
		}
	}

	for (int gamma = 0; gamma < Gamma; gamma++) {
		if (BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {

			//Se tem em excesso, remover dos padrões que só tem o tipo gamma
			for (int i = 0; i < n_genes; i++) {
				if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
					int indice_real = solu.ind[i] - P;

					//só gera formas do tipo gamma?
					bool gera_gamma = true;
					for (int tam = 0; tam < Gamma; tam++) {
						if (tam != gamma && CutPatterns[indice_real].tamanhos[tam] > 0) {
							gera_gamma = false;
							break;
						}
					}
					if (gera_gamma && CutPatterns[indice_real].tamanhos[gamma] == 0)
						gera_gamma = false;

					if (gera_gamma) {
						int remover = min(solu.n_vezes[i], (int)ceil((double)(BarrasGeradas[gamma] - FormasQueDevemSerGeradas[gamma])
							/ CutPatterns[indice_real].tamanhos[gamma]));

						solu.n_vezes[i] -= remover;
						BarrasGeradas[gamma] -= remover * CutPatterns[indice_real].tamanhos[gamma];
						EstoqueUsado[CutPatterns[indice_real].index_barra] -= remover;


					}
					if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
						break;
				}

			}


			//Se continuar, remova dos traspasses
			if (BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {
				for (int i = 0; i < n_genes; i++){
					if (solu.ind[i] >= P + H && solu.ind[i] <= P - 1 + H + O) {
						int indice_real = solu.ind[i] - P - H;

						if (SplPatterns[indice_real].barra_gerada == gamma && solu.n_vezes[i] > 0) {
							int remover = min(solu.n_vezes[i], BarrasGeradas[gamma] - FormasQueDevemSerGeradas[gamma]);
							solu.n_vezes[i] -= remover;
							BarrasGeradas[gamma] -= remover;

							for (int v = 0; v < V; v++) {
								EstoqueUsado[W + v] -= remover * SplPatterns[indice_real].tamanhos[v];
							}
						}
						if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
							break;
					}
				}
			}
		}



		if (BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma]) {

			//Se ta em falta, adiciona
			for (int i = 0; i < n_genes; i++) {
				if (solu.ind[i] >= P && solu.ind[i] <= P - 1 + H) {
					int indice_real = solu.ind[i] - P;
					//Se gera só gamma
					bool gera_gamma = true;
					for (int tam = 0; tam < Gamma; tam++) {
						if (tam != gamma && CutPatterns[indice_real].tamanhos[tam] > 0) {
							gera_gamma = false;
							break;
						}
					}
					if (gera_gamma && CutPatterns[indice_real].tamanhos[gamma] == 0)
						gera_gamma = false;


					//se o padrão cobre uma barra de tamanho gamma
					if (gera_gamma) {
						//enquanto ainda é menor e o padrão não interfere no estoque da barra
						int adicionar = min((int)(floor((double)(FormasQueDevemSerGeradas[gamma] - BarrasGeradas[gamma]) / CutPatterns[indice_real].tamanhos[gamma])),
							e[CutPatterns[indice_real].index_barra] - EstoqueUsado[CutPatterns[indice_real].index_barra]);

						EstoqueUsado[CutPatterns[indice_real].index_barra] += adicionar;
						solu.n_vezes[i] += adicionar;

						BarrasGeradas[gamma] += adicionar * CutPatterns[indice_real].tamanhos[gamma];

						if (BarrasGeradas[gamma] == FormasQueDevemSerGeradas[gamma])
							break;
					}
				}
			}
			//Se ainda não supriu adiciona splicing feito doido
			if (BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma]) {
				for (int i = 0; i < n_genes; i++) {
					if (solu.ind[i] >= P + H && solu.ind[i] <= P - 1 + H + O) {
						int indice_real = solu.ind[i] - P - H;
						if (SplPatterns[indice_real].barra_gerada == gamma) {
							bool pode_ser_adicionado = true;
							while (pode_ser_adicionado  && BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma]) {

								for (int v = 0; v < V; v++) {
									if (EstoqueUsado[W + v] + SplPatterns[indice_real].tamanhos[v] > e[W + v]) {
										pode_ser_adicionado = false;
										break;
									}
								}
								if (pode_ser_adicionado) {
									solu.n_vezes[i]++;
									BarrasGeradas[gamma]++;
									for (int v = 0; v < V; v++)
										EstoqueUsado[W + v] += SplPatterns[indice_real].tamanhos[v];
								}
							}
							if (BarrasGeradas[gamma] == FormasQueDevemSerGeradas[gamma])
								break;
						}
					}
				}
			}

		}
	}


	/*if (!viavel(solu)) {
		cout << "Erro!" << endl;
		system("pause");
		exit(1);
	}*/

	//Remover os elementos com n_vezes iguais a zero
	vector<int> serao_removidos(0, 0);
	for (int i = 0; i < n_genes; i++) {
		if (solu.n_vezes[i] == 0) {
			serao_removidos.push_back(i);
		}
	}
	for (int i = serao_removidos.size() - 1; i >= 0; i--) {
		solu.ind.erase(solu.ind.begin() + serao_removidos[i]);
		solu.n_vezes.erase(solu.n_vezes.begin() + serao_removidos[i]);
	}
}

GA_Novo::individuo GA_Novo::cruzamento_unico(vector<individuo> Populacao)
{
	uniform_int_distribution<int> distribuicao(0, Populacao.size() - 1);

	int i = distribuicao(generator),
		j = distribuicao(generator);
	while (i == j)
		j = distribuicao(generator);


	return cruzar(Populacao[i], Populacao[j]);
}

void GA_Novo::mutar(individuo &solu){
	int n_genes = solu.ind.size();
	
	uniform_int_distribution<int> 
		uniformeP(1, P - 1 + H + O),
		uniforme_mutar(0, n_genes - 1);

	int mutado;
	do
	{
		mutado = uniforme_mutar(generator);
	} while (solu.n_vezes[mutado] == 0);
	solu.n_vezes[mutado] --;



	int mutante;
	bool ja_esta;
	do {
		mutante = uniformeP(generator);
		ja_esta = (find(solu.ind.begin(), solu.ind.end(), mutante) != solu.ind.end());
	} while (ja_esta);

	solu.ind.push_back(mutante);
	solu.n_vezes.push_back(1);
	corrigir(solu);
}

void GA_Novo::Restart(vector<individuo>& Populacao)
{
	//Orderna populacao por fitness
	std::sort(Populacao.begin(), Populacao.end(), [](individuo i1, individuo i2) {return i1.fitness < i2.fitness; });
	//Salva populacao antiga
	vector<individuo> Auxiliar = Populacao;
	//Limpar populacao para adicionar os individuos selecionados
	Populacao.clear();

	//Elistimo
	for(int i = 0; i < 10; i++)
		Populacao.push_back(Auxiliar[i]);
	
	for (int i = 0; i < max(TamanhoDaPopulacao*10, 10*(P + O + H)); i++)
	{
		individuo solucao = GerarSoluGRASP();

		if (!viavel(solucao))
			corrigir(solucao);

		if (viavel(solucao)) {
			solucao.fitness = fitness(solucao);
			Populacao.push_back(solucao);
		}

	}
	selecao(Populacao);
}


void GA_Novo::ImprimirArquivo(individuo solu, double time)
{
	ofstream resultados("resultados_GA.txt", fstream::app);

	resultados << endl << nome_instancia << "," << solu.fitness << "," << time << endl;

	resultados.close();
}