#include "Modelo_Cplex.h"
using namespace std;


int main() {
	
	Modelo_Cplex Modelo("arquivoteste2");

	Modelo.MontarModelo();
	Modelo.resolver_inteira();
	Modelo.ImprimirSolucao();
	Modelo.ImprimirGantt();

	getchar();
	return 0;
}

