#include "Modelo_Cplex.h"
using namespace std;


int main() {
	
	Modelo_Cplex Modelo("problema");

	Modelo.MontarModelo();
	Modelo.resolver_linear();
	getchar();
	return 0;
}

