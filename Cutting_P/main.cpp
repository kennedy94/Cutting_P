#include "Modelo_Cplex.h"
#include "Heuristica.h"
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

	//Modelo.ImprimirGantt();
	//Modelo.PlotarBarras();
	
	Heuristica GA(inst);
	GA.funcaoteste();
	return 0;
}

