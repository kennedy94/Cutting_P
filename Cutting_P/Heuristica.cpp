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
			//filho.ind[i] = mae.ind[i];
			filho.n_vezes[i] = mae.n_vezes[i];
		}
	if (!viavel(filho))
		corrigir(filho);
	
	filho.fitness = fitness(filho);
	filhos.push_back(filho);


	//filho2
	filho = pai;
	for (int i = 0; i < filho.ind.size(); i++)
		if (i % 2 != 0 ) {
			//filho.ind[i] = mae.ind[i];
			filho.n_vezes[i] = mae.n_vezes[i];
		}
	if (!viavel(filho))
		corrigir(filho);


	filho.fitness = fitness(filho);
	filhos.push_back(filho);

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
	}


	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++)
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
			if (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k])
				return false;
			
	for (int i = 0; i < P - 1; i++)
	{
		for (int gamma = 0; gamma < Gamma; gamma++)
			if (PackPatterns[solu.ind[i]].maximal(L[gamma], TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos))
				FormasQueDevemSerGeradas[gamma] += TipoVigas[PackPatterns[solu.ind[i]].tipo].n_barras*solu.n_vezes[i];
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
			//cout << endl << endl;
			//ImprimirVetorSolu(solu);
			return false;
		}
	}



	vector<int> BarrasGeradasPelaSolu(Gamma, 0);
	//contar as barras com tamanhos de forma geradas
	for (int gamma = 0; gamma < Gamma; gamma++) {
		for (int i = P - 1; i < P - 1 + H + O; i++) {
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












void Heuristica::corrigir(individuo &solu)
{
	
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

				while (DemandasAuxiliares[c].demandas[k] < TipoVigas[c].demandas[k]) {
					solu.n_vezes[primeiro_padrao]++;
					for (int k2 = 0; k2 < TipoVigas[c].n_comprimentos; k2++)
						DemandasAuxiliares[c].demandas[k2] += PackPatterns[solu.ind[primeiro_padrao]].tamanhos[k2];
				}
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
							int(floor(retirar / SplPatterns[solu.ind[i]].tamanhos[w - W])));
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
					while (solu.n_vezes[i] > 0 && BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {
						solu.n_vezes[i]--;
						//Como estou assumindo que o padrão selecionado só gera barra gamma então só basta atualizar o valor de gamma
						BarrasGeradas[gamma] -= CutPatterns[solu.ind[i]].tamanhos[gamma];
						EstoqueUsado[CutPatterns[solu.ind[i]].index_barra]--;
					}
					if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
						break;
				}
			}


			//Se continuar, remova dos traspasses
			if (BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {
				for (int i = P - 1 + H; i < P - 1 + H + O; i++) {
					if (SplPatterns[solu.ind[i]].barra_gerada == gamma) {
						while (solu.n_vezes[i] > 0 && BarrasGeradas[gamma] > FormasQueDevemSerGeradas[gamma]) {
							solu.n_vezes[i]--;
							BarrasGeradas[gamma]--;

							for (int v = 0; v < V; v++)
								EstoqueUsado[W + v] -= SplPatterns[solu.ind[i]].tamanhos[v];

						}
						if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
							break;
					}
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
					while (BarrasGeradas[gamma] < FormasQueDevemSerGeradas[gamma] 
						&& EstoqueUsado[CutPatterns[solu.ind[i]].index_barra] < e[CutPatterns[solu.ind[i]].index_barra]) {
						
						EstoqueUsado[CutPatterns[solu.ind[i]].index_barra]++;
						solu.n_vezes[i]++;

						BarrasGeradas[gamma] += CutPatterns[solu.ind[i]].tamanhos[gamma];
					}
					if (BarrasGeradas[gamma] <= FormasQueDevemSerGeradas[gamma])
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


	vector<int> BarrasGeradasPelaSolu(Gamma, 0);
	//contar as barras com tamanhos de forma geradas
	for (int gamma = 0; gamma < Gamma; gamma++) {
		for (int i = P - 1; i < P - 1 + H + O; i++) {
			if (i < P - 1 + H)
				BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
			else
				if (SplPatterns[solu.ind[i]].barra_gerada == gamma)
					BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i];
		}
	}
	
	for (int gamma = 0; gamma < Gamma; gamma++)
		if (BarrasGeradasPelaSolu[gamma] != FormasQueDevemSerGeradas[gamma])
			cout << "CUIDADOOO" << endl;


	//Corrigir por necessidade de barras
	if (!viavel(solu)) {
		cout << "help" << endl;
		viavel(solu);
	}
	else
		cout << "Corrigiu!" << endl;
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
				int n_vezes_uso_padrao_atual = ceil(qtde_necessaria[k_aux] / PackPatterns[solu1[i]].tamanhos[k_aux]);

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
				//cout << endl << "gamma:" << gamma << " - " <<
				//	EstoqueDisponivel[CutPatterns[solu1[i]].index_barra] << " , " << necessario << endl << endl;
				
				int adicionado = floor(min(EstoqueDisponivel[CutPatterns[solu1[i]].index_barra], necessario) / CutPatterns[solu1[i]].tamanhos[gamma]);
				FormasGeradas[gamma] += CutPatterns[solu1[i]].tamanhos[gamma] * adicionado;
				solu2[i] += adicionado;
				EstoqueDisponivel[CutPatterns[solu1[i]].index_barra] -= adicionado;
				necessario -= CutPatterns[solu1[i]].tamanhos[gamma] * adicionado;
			}
		}
		//cout << "nec:" << necessario << endl;
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
						if(adicionado > floor(EstoqueDisponivel[W + v] / SplPatterns[solu1[i]].tamanhos[v]))
							adicionado = floor(EstoqueDisponivel[W + v] / SplPatterns[solu1[i]].tamanhos[v]);
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
	srand(time(NULL));

	double MenorFo = double(INT_MAX);
	

	vector<individuo> Populacao;


	cout << "Gerando populacao inicial\n";

	auto start = chrono::system_clock::now();

	for (int i = 0; i < 100; i++)
	{
		individuo solucao = GerarSoluGRASP();

		if (!viavel(solucao)) {
			cout << "Erro! Solução gulosa gerada é inviável!" << endl;
			system("pause");
			exit(0);
		}

		solucao.fitness = fitness(solucao);

		Populacao.push_back(solucao);

		if ( (i + 1) % 100 == 0) {
			auto end = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
			cout << "\t" << setprecision(2) << elapsed_seconds.count() << "s" << endl;
		}
	}
	
	selecao(Populacao);

	list<individuo> Offspring = cruzamento(Populacao);
	Populacao.insert(Populacao.end(), Offspring.begin(), Offspring.end());
	selecao(Populacao);


	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "\t Tempo Total " << setprecision(2) << elapsed_seconds.count() << "s" << endl;

	cout << endl << "Populacao inicial gerada com sucesso" << endl;
	cout << viavel(Populacao[0]) << endl;
	cout << "FO = " << Populacao[0].fitness << endl;
}