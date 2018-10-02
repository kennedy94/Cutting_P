#include "Heuristica.h"

// Function to return k'th smallest element in a given array
inline int kMenor(const vector<int> v, int k) {

	// initialize original index locations
	vector<int> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] < v[i2]; });

	return idx[k - 1];
}


inline int kMaior(const vector<int> v, int k) {

	// initialize original index locations
	vector<int> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2]; });

	return idx[k - 1];
}




double Heuristica::fitness(individuo solu)
{
	vector<int> UtilizacaoFormas(M, 0);
	//vector<vector<int>> PadroesAssociados(M);

	int Makespan;
	double 
		Sobra1 = 0.0,
		Sobra2 = 0.0,
		Sobra3 = 0.0;

	/*for (int m = 0; m < M; m++) 
		PadroesAssociados[m].resize(T);*/
	
	//Passa nos genes dependendes dos padroes de empacotamento
	//é -1 porque não tem indíce 0 na heuristica
	for (int i = 0; i < P-1; i++)
	{
		for (int vezes = 0; vezes < solu.n_vezes[i]; vezes++) {

			//forma atual recebe argmin da utilização que acomoda o padrão
			int k_esimo = 0;
			int FormaAtual;
			do
			{
				k_esimo++;
				FormaAtual = kMenor(UtilizacaoFormas, k_esimo);
			} while (!PackPatterns[solu.ind[i]].maximal(FORMAS[FormaAtual],TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos));
			
			//tempo atual é o valor mínimo da utilização
			int TempoAtual = UtilizacaoFormas[kMenor(UtilizacaoFormas, k_esimo)];

			//PadroesAssociados[FormaAtual][TempoAtual] = solu.ind[i];

			UtilizacaoFormas[FormaAtual] += TipoVigas[PackPatterns[solu.ind[i]].tipo].tempo_cura;
		}
	}

	Makespan = *max_element(UtilizacaoFormas.begin(), UtilizacaoFormas.end());

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

	//filho1
	individuo filho = pai;
	for(int i = 0; i < filho.ind.size(); i++)
		if (i % 2 == 0) {
			filho.ind[i] = mae.ind[i];
			filho.n_vezes[i] = mae.n_vezes[i];
		}
	//if (!viavel(filho))
		corrigir(filho);
	filho.fitness = fitness(filho);
	filhos.push_back(filho);


	//filho2
	filho = pai;
	for (int i = 0; i < filho.ind.size(); i++)
		if (i % 2 != 0 ) {
			filho.ind[i] = mae.ind[i];
			filho.n_vezes[i] = mae.n_vezes[i];
		}
	//if (!viavel(filho))
		corrigir(filho);
	filho.fitness = fitness(filho);
	filhos.push_back(filho);


	//filho3
	/*filho = pai;
	for (int i = ceil((P - 1)/2); i < P - 1; i++){
		filho.ind[i] = mae.ind[i];
		filho.n_vezes[i] = mae.n_vezes[i];
	}
	for (int i = ceil((P - 1 + H)/2); i < P - 1 + H; i++) {
		filho.ind[i] = mae.ind[i];
		filho.n_vezes[i] = mae.n_vezes[i];
	}

	for (int i = ceil((P - 1 + H + O) / 2); i < P - 1 + H + O; i++) {
		filho.ind[i] = mae.ind[i];
		filho.n_vezes[i] = mae.n_vezes[i];
	}
	if (!viavel(filho))
		corrigir(filho);
	filho.fitness = fitness(filho);
	filhos.push_back(filho);


	//filho4
	filho = pai;
	for (int i = 0; i < ceil((P - 1) / 2); i++) {
		filho.ind[i] = mae.ind[i];
		filho.n_vezes[i] = mae.n_vezes[i];
	}
	for (int i = P - 1; i < ceil((P - 1 + H) / 2); i++) {
		filho.ind[i] = mae.ind[i];
		filho.n_vezes[i] = mae.n_vezes[i];
	}

	for (int i = P - 1 + H; i < ceil((P - 1 + H + O) / 2); i++) {
		filho.ind[i] = mae.ind[i];
		filho.n_vezes[i] = mae.n_vezes[i];
	}
	if (!viavel(filho))
		corrigir(filho);
	filho.fitness = fitness(filho);
	filhos.push_back(filho);*/

	return filhos;
}





list<individuo> Heuristica::cruzamento(vector<individuo> Popu)
{
	list<individuo> retorno;
	for (int i = 0; i < Popu.size(); i++) {
		for (int j = 0; j < Popu.size(); j++) {
			if (i != j) {
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

		for (int gamma = 0; gamma < Gamma; gamma++)
			if (PackPatterns[solu.ind[i]].maximal(L[gamma], TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos))
				FormasQueDevemSerGeradas[gamma] += TipoVigas[PackPatterns[solu.ind[i]].tipo].n_barras*solu.n_vezes[i];
	}


	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++)
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
			if (DemandasAuxiliares[c].comprimentos[k] < TipoVigas[c].comprimentos[k])
				return false;
			
	//Viabilidade por estoque

	vector<int> EstoqueUsado(W + V, 0);
	
	for (int i = P - 1; i < P - 1 + H + O; i++) {
		if (i < P - 1 + H)
			EstoqueUsado[CutPatterns[solu.ind[i]].index_barra] += solu.n_vezes[i];
		else {
			for (int v = 0; v < V; v++)
				EstoqueUsado[W + v] += SplPatterns[solu.ind[i]].tamanhos[v];
		}
	}

	for(int w = 0; w < W+V; w++)
		if (EstoqueUsado[w] > e[w])
			return false;



	vector<int> BarrasGeradasPelaSolu(Gamma, 0);
	//contar as barras com tamanhos de forma geradas
	for (int i = P - 1; i < P - 1 + H + O; i++) {
		for (int gamma = 0; gamma < Gamma; gamma++) {
			if (i < P - 1 + H)
				BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
			else
				if (SplPatterns[solu.ind[i]].barra_gerada == gamma)
					BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i];
		}
	}

	//Se a quantidade da solu for diferente da que deveria ser é inviável
	for (int gamma = 0; gamma < Gamma; gamma++)
		if (BarrasGeradasPelaSolu[gamma] != FormasQueDevemSerGeradas[gamma])
			return false;

	return true;
}












void Heuristica::corrigir(individuo & solu)
{
	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;

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

		//Contar o número de formas de tamanho L[gamma] utilizadas
		for (int gamma = 0; gamma < Gamma; gamma++)
			if (PackPatterns[solu.ind[i]].maximal(L[gamma], TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos))
				FormasQueDevemSerGeradas[gamma] += TipoVigas[PackPatterns[solu.ind[i]].tipo].n_barras*solu.n_vezes[i];
	}





	bool corrigiu = false;
	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++)
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
			if (DemandasAuxiliares[c].comprimentos[k] < TipoVigas[c].comprimentos[k]) {
				int primeiro_padrao;
				for (int i = 0; i < P - 1; i++) {
					if (PackPatterns[solu.ind[i]].tipo == c && PackPatterns[solu.ind[i]].tamanhos[k] > 0 && solu.n_vezes[i]) {
						primeiro_padrao = i;
						break;
					}
				}

				while (DemandasAuxiliares[c].comprimentos[k] < TipoVigas[c].comprimentos[k]) {
					solu.n_vezes[primeiro_padrao]++;

					for (int k2 = 0; k2 < TipoVigas[c].n_comprimentos; k2++)
						DemandasAuxiliares[c].comprimentos[k2] += PackPatterns[solu.ind[primeiro_padrao]].tamanhos[k2];
				}
				corrigiu = true;
			}

	if (corrigiu)
		cout << "Corrigiu por demanda de viga" << endl;

			
	//Viabilidade por estoque

	vector<int> EstoqueUsado(W + V, 0);

	for (int w = 0; w < W + V; w++) {
		for (int i = P - 1; i < P - 1 + H + O; i++) {
			if (i < P - 1 + H && CutPatterns[solu.ind[i]].index_barra == w)
				EstoqueUsado[w] += solu.n_vezes[i];
			else
				EstoqueUsado[w] += SplPatterns[solu.ind[i]].tamanhos[w];
		}
		if (EstoqueUsado[w] > e[w]) {
			//Vou tratar isso nos operadores, então o ideal é não gerar solução desse tipo
			cerr << "Erro solução inviável por estoque!" << endl;
			exit(0);
		}
	}

	vector<int> BarrasGeradasPelaSolu(Gamma, 0);
	//contar as barras com tamanhos de forma geradas
	for (int i = P - 1; i < P - 1 + H + O; i++) {
		for (int gamma = 0; gamma < Gamma; gamma++) {
			if (i < P - 1 + H)
				BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
			else
				if (SplPatterns[solu.ind[i]].barra_gerada == gamma)
					BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i];
		}
	}

	//Se a quantidade da solu for diferente da que deveria ser é inviável
	for (int gamma = 0; gamma < Gamma; gamma++) {
		list<int> proibidos;
		while (BarrasGeradasPelaSolu[gamma] != FormasQueDevemSerGeradas[gamma]) {
			int primeiro_padrao = -1;

			//preferencia para os padrões que já estão na solução
			for (int i = P - 1; i < P - 1 + H; i++) {
				if (CutPatterns[solu.ind[i]].tamanhos[gamma] > 0 && solu.n_vezes[i] > 0) {
					if (proibidos.empty()) {
						primeiro_padrao = i;
						break;
					}
					else {
						if (proibidos.end() == find(proibidos.begin(), proibidos.end(), i)) {
							primeiro_padrao = i;
							break;
						}
					}
				}
			}
			
			//se não houver padrões que cubram a barra gamma na solução atual pegar o primeiro padrão que cubra
			if (primeiro_padrao == -1) {
				for (int i = P - 1; i < P - 1 + H; i++) {
					if (CutPatterns[solu.ind[i]].tamanhos[gamma] > 0) {
						if (proibidos.empty()) {
							primeiro_padrao = i;
							break;
						}
						else {
							if (proibidos.end() == find(proibidos.begin(), proibidos.end(), i)) {
								primeiro_padrao = i;
								break;
							}
							else
							{
								continue;
							}
							
						}
					}
				}
				if (primeiro_padrao == -1)//não tem estoque, adiciona traspasse
				{
					for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
						if (SplPatterns[solu.ind[i]].barra_gerada == gamma) {
							list<int> usados;
							for (int tam = 0; tam < SplPatterns[solu.ind[i]].tamanhos.size(); tam++) {
								SplPatterns[solu.ind[i]].tamanhos[tam] > 0;
								usados.push_back(W+tam);
							}

							while (BarrasGeradasPelaSolu[gamma] < FormasQueDevemSerGeradas[gamma] && estoque_ok(EstoqueUsado,usados)) {
								BarrasGeradasPelaSolu[gamma]++;
								solu.n_vezes[i]++;
								for (int tam = 0; tam < SplPatterns[solu.ind[i]].tamanhos.size(); tam++)
									EstoqueUsado[W + tam] += SplPatterns[solu.ind[i]].tamanhos[tam];

							}
						}
					}
				}
			}

			while (primeiro_padrao != -1 && BarrasGeradasPelaSolu[gamma] < FormasQueDevemSerGeradas[gamma] &&
				EstoqueUsado[CutPatterns[solu.ind[primeiro_padrao]].index_barra] < e[CutPatterns[solu.ind[primeiro_padrao]].index_barra])//se é inviável e o estoque não é excedido
			{
				//adicionar restrições de estoque aqui
				solu.n_vezes[primeiro_padrao]++;
				BarrasGeradasPelaSolu[gamma] += CutPatterns[solu.ind[primeiro_padrao]].tamanhos[gamma];
				EstoqueUsado[CutPatterns[solu.ind[primeiro_padrao]].index_barra]++;
					
			}

			while (BarrasGeradasPelaSolu[gamma] > FormasQueDevemSerGeradas[gamma] && solu.n_vezes[primeiro_padrao] > 0) {
				solu.n_vezes[primeiro_padrao]--;
				BarrasGeradasPelaSolu[gamma] -= CutPatterns[solu.ind[primeiro_padrao]].tamanhos[gamma];
			}
			
			proibidos.push_back(primeiro_padrao);
		}
	}
}





bool Heuristica::estoque_ok(vector<int> EstoqueUsado, list<int> usados) {
	for (auto w: usados){
		if (EstoqueUsado[w] >= e[w])
			return false;
	}
	return true;
}






void Heuristica::ImprimirVetorSolu(individuo solu)
{
	for (int i = P - 1; i < P - 1 + H + O; i++) {
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
	/*random_shuffle(solu1.begin(), solu1.begin() + P - 1);
	random_shuffle(solu1.begin() + P - 1, solu1.begin() + P - 1 + H);
	random_shuffle(solu1.begin() + P - 1 + H, solu1.end());*/


	//prenchendo os padrões de empacotamento até a demanda ser atendida
	vector<Tipo_Viga> DemandasAuxiliares = TipoVigas;

	//Iniciar as demandas atuais em 0
	for (auto &elemento : DemandasAuxiliares)
		for (auto &demand : elemento.demandas)
			demand = 0;

	for (int i = 0; i < P - 1; i++) {
		int tipo_atual = PackPatterns[solu1[i]].tipo;
		
		bool demanda_total_atendida = true;
		for (int c = 0; c < C; c++)
			for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
				if (DemandasAuxiliares[tipo_atual].demandas[k] < TipoVigas[tipo_atual].demandas[k])
					demanda_total_atendida = false;
		if (demanda_total_atendida) 
			break;


		bool demanda_tipo_atendida = true;
		vector<int> qtde_necessaria(TipoVigas[tipo_atual].n_comprimentos, 0);
		for (int k = 0; k < TipoVigas[tipo_atual].n_comprimentos; k++) {
			if (DemandasAuxiliares[tipo_atual].demandas[k] < TipoVigas[tipo_atual].demandas[k]) {
				demanda_tipo_atendida = false;
				qtde_necessaria[k] = TipoVigas[tipo_atual].demandas[k] - DemandasAuxiliares[tipo_atual].demandas[k];
			}
		}
		demanda_tipo_atendida;
		if (!demanda_tipo_atendida) {
			int maior_necessaria = -1;
			//teste se o padrão cobre o tamanho necessario
			//pega o tamanho mais necessário que o padrão cobre
			for (int kmaior = 1; kmaior <= TipoVigas[tipo_atual].n_comprimentos; kmaior++) {
				int tamanho_escolhido = kMaior(qtde_necessaria, kmaior);
				if (PackPatterns[solu1[i]].tamanhos[tamanho_escolhido] > 0)
					maior_necessaria = tamanho_escolhido;
			}

			if (maior_necessaria == -1) continue;
			int n_vezes_uso_padrao_atual = ceil(qtde_necessaria[maior_necessaria] / PackPatterns[solu1[i]].tamanhos[maior_necessaria]);
			//atualizar as demandas auxiliares
			solu2[i] = n_vezes_uso_padrao_atual;

			for (int k = 0; k < TipoVigas[tipo_atual].n_comprimentos; k++) {
				DemandasAuxiliares[tipo_atual].demandas[k] +=
					n_vezes_uso_padrao_atual*PackPatterns[solu1[i]].tamanhos[k];
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
	cout << "\n\n\n";
	for (int gamma = 0; gamma < Gamma; gamma++) {
		int necessario = FormasQueDevemSerGeradas[gamma];
		for (int i = P - 1; i < P - 1 + H; i++) {
			if (CutPatterns[solu1[i]].tamanhos[gamma] > 0) {
				int adicionado = min(EstoqueDisponivel[CutPatterns[solu1[i]].index_barra], necessario);
				FormasGeradas[gamma] += CutPatterns[solu1[i]].tamanhos[gamma]* adicionado;
				solu2[i] += adicionado;
				EstoqueDisponivel[CutPatterns[solu1[i]].index_barra] -= adicionado;
				necessario -= CutPatterns[solu1[i]].tamanhos[gamma] * adicionado;
			}
		}
		
	}
	//cout << "\n o/ " << viavel(solucao) << endl;
	individuo solucao;
	solucao.ind = solu1;
	solucao.n_vezes = solu2;
	cout << "\n\n\nGERADA\n";
	ImprimirVetorSolu(solucao);


	viavel(solucao);

	return solucao;
}







individuo Heuristica::GerarSoluAleatoria()
{

	vector<int>	solu1(P - 1 + H + O),	//índices dos padrões
				solu2(P - 1 + H + O);	//quantidade dos padrões

								//Inicializando os índices
	iota(solu1.begin(), solu1.begin() + P - 1, 1);
	iota(solu1.begin() + P - 1, solu1.begin() + P - 1 + H, 0);
	iota(solu1.begin() + P - 1 + H, solu1.end(), 0);

	//Definição de sementes para a geração das quantidades
	mt19937 rng_int, rng_real;
	rng_int.seed(random_device()());
	rng_real.seed(random_device()());

	int SomaPack = 0,
		SomaCut = 0,
		SomaSpl = 0;

	for (auto i = solu2.begin(); i != solu2.begin() + P - 1; i++) {
		uniform_real_distribution<double> u_real(0.0, 1.0);

		if (u_real(rng_real) < 0.8) {
			*i = 0;
		}
		else {
			uniform_int_distribution<int> u_inteira(1, ceil(T / M));
			*i = u_inteira(rng_int);
			SomaPack += *i;
		}
	}

	int w = P - 1;
	for (auto i = solu2.begin() + P - 1; i != solu2.begin() + P - 1 + H; i++) {
		uniform_real_distribution<double> u_real(0.0, 1.0);

		if (u_real(rng_real) < 0.5) {
			*i = 0;
		}
		else {
			uniform_int_distribution<int> u_inteira(0, min(e[CutPatterns[solu1[w]].index_barra], int(SomaPack / H)));
			*i = u_inteira(rng_int);
		}
		w++;
	}

	for (auto i = solu2.begin() + P - 1; i != solu2.begin() + P - 1 + H; i++) {
		uniform_int_distribution<int> u_inteira(0, int(SomaPack / O));
		*i = u_inteira(rng_int);
	}


	//cout << ":)" << endl;
	individuo solucao;
	solucao.ind = solu1;
	solucao.n_vezes = solu2;

	return solucao;
}



void Heuristica::selecao(vector<individuo> &Popu)	//Selecao por elitismo
{
	//Ordenando população por fitness
	sort(Popu.begin(), Popu.end(), [](individuo i1, individuo i2) {return i1.fitness < i2.fitness; });

	TamanhoDaPopulacao = 10;
	//Auxiliar para guardar a populacao atual
	vector<individuo> Auxiliar = Popu;
	Popu.clear();
	//Populacao atual recebe só os TamanhoDaPopulacao melhores elementos
	for (int i = 0; i < TamanhoDaPopulacao; i++)
		Popu.push_back(Auxiliar[i]);
}





void Heuristica::funcaoteste() {

	//GerarSoluGRASP();

	double MenorFo = double(INT_MAX);

	vector<individuo> Populacao;


	cout << "Gerando populacao inicial\n";

	auto start = chrono::system_clock::now();

	for (int i = 0; i < 200; i++)
	{
		individuo solucao = GerarSoluGRASP();

		//corrigir(solucao);
			
		solucao.fitness = fitness(solucao);

		if (MenorFo > solucao.fitness)
			MenorFo = solucao.fitness;

		//if(!viavel(solucao))
		//	cout << "0 ";

		

		Populacao.push_back(solucao);
		if (i % 100 == 0) {
			auto end = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
			cout << "\t" << elapsed_seconds.count() << "s" << endl;
		}
	}
	selecao(Populacao);

	list<individuo> Offspring = cruzamento(Populacao);
	Populacao.insert(Populacao.end(), Offspring.begin(), Offspring.end());



	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "\t Tempo Total " << elapsed_seconds.count() << "s" << endl;
	
	selecao(Populacao);

	cout << endl << "Populacao inicial gerada com sucesso" << endl;

	cout << "FO = " << MenorFo << endl;
}