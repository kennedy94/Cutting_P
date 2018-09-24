#include "Heuristica.h"

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

	Makespan = *max(UtilizacaoFormas.begin(), UtilizacaoFormas.end());

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


	return filhos;
}


// Function to return k'th smallest element in a given array
int kMenor(const vector<int> v, int k) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx[k-1];
}

bool Heuristica::viavel(individuo solu)
{
	vector<int> FormasQueDevemSerGeradas(Gamma, 0);
	vector<int> UtilizacaoFormas(M, 0);

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


		for (int vezes = 0; vezes < solu.n_vezes[i]; vezes++) {
			//forma atual recebe argmin da utilização que acomoda o padrão
			int k_esimo = 0;
			int FormaAtual;
			do
			{
				k_esimo++;
				FormaAtual = kMenor(UtilizacaoFormas, k_esimo);
			} while (!PackPatterns[solu.ind[i]].maximal(FORMAS[FormaAtual], TipoVigas[PackPatterns[solu.ind[i]].tipo].comprimentos));

			//tempo atual é o valor mínimo da utilização
			int TempoAtual = UtilizacaoFormas[kMenor(UtilizacaoFormas, k_esimo)];

			//Contar o número de formas de tamanho L[gamma] utilizadas
			for (int gamma = 0; gamma < Gamma; gamma++)
				if (FORMAS[FormaAtual] == L[gamma])
					FormasQueDevemSerGeradas[gamma]++;

			UtilizacaoFormas[FormaAtual] += TipoVigas[PackPatterns[solu.ind[i]].tipo].tempo_cura;
		}
	}

	//Se não obedece a demanda de pelo menos uma viga é inviável por Demanda de Viga
	for (int c = 0; c < C; c++)
		for (int k = 0; k < TipoVigas[c].n_comprimentos; k++)
			if (DemandasAuxiliares[c].comprimentos[k] < TipoVigas[c].comprimentos[k])
				return false;

	vector<int> BarrasGeradasPelaSolu(Gamma, 0);
	//contar as barras com tamanhos de forma geradas
	for (int i = P - 1; i < P - 1 + H + O; i++) {
		for (int gamma = 0; gamma < Gamma; gamma++){
			if (i < P - 1 + H)
				BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i] * CutPatterns[solu.ind[i]].tamanhos[gamma];
			else
				if (CutPatterns[solu.ind[i]].index_barra == gamma)
					BarrasGeradasPelaSolu[gamma] += solu.n_vezes[i];
		}	
	}

	//Se a quantidade da solu for diferente da que deveria ser é inviável
	for (int gamma = 0; gamma < Gamma; gamma++)
		if (BarrasGeradasPelaSolu[gamma] != FormasQueDevemSerGeradas[gamma])
			return false;

	//Viabilidade por estoque

	vector<int> EstoqueUsado(W + V, 0);

	for (int w = 0; w < W + V; w++) {
		for (int i = P - 1; i < P - 1 + H + O; i++) {
			if (i < P - 1 + H && CutPatterns[solu.ind[i]].index_barra == w)
				EstoqueUsado[w] += solu.n_vezes[i];
			else
				EstoqueUsado[w] += SplPatterns[solu.ind[i]].tamanhos[w];
		}
		if (EstoqueUsado[w] > e[w])
			return false;

	}

	return true;
}

void Heuristica::funcaoteste(){

	vector<int> Teste = { 50, 20, 60 };

	cout << kMenor(Teste, 1) << endl;

}
