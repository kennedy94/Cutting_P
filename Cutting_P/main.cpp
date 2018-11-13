#include "Modelo_Cplex.h"
#include "Heuristica.h"
#include "GA_Novo.h"
#include "Relaxacao_Linear.h"
using namespace std;


int main(int argc, char *argv[]) {

	char *inst;
	if (argc < 2)
		inst = "problema.txt";
	else {
		if (argc < 3)
			inst = argv[1];
		else {
			cout << "Argumentos demais" << endl;
			exit(0);
		}
	}

	/*Modelo_Cplex Modelo(inst);
	Modelo.MontarModelo();
	Modelo.resolver_inteira();
	Modelo.ImprimirSolucaoArquivo();
	Modelo.ImprimirSolucao();*/

	/*Relaxacao_Linear Linear(inst);
	Linear.MontarModelo();
	Linear.resolver_linear();
	Linear.ImprimirSolucaoArquivo();*/

	
	//Modelo.ImprimirGantt();
	//Modelo.PlotarBarras();

	/*for (int i = 0; i < 1; i++){
	Heuristica GA(inst);
	GA.funcaoteste();
	}*/

	for (int i = 0; i < 1; i++) {
		GA_Novo GA(inst);
		GA.funcao_teste();
	}

	return 0;
}

