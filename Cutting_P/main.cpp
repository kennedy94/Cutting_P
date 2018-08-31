#include "Modelo_Cplex.h"
using namespace std;


int main() {
	
	Modelo_Cplex Modelo("cwp1");

	Modelo.MontarModelo();
	Modelo.resolver_inteira();
	Modelo.ImprimirSolucao();
	Modelo.ImprimirGantt();
	Modelo.PlotarBarras();

	getchar();
	return 0;
}

