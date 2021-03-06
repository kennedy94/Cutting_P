%Entrada: 
%   C = N�mero de tipos de vigas
%   M = N�mero de formas
%   arq = Nome do arquivo a ser gerado. Ex: 'inst1.txt'
%   misturado = 1 se tempo de cura s�o misturados, 0 para usar o padr�o
%       Por padr�o o tipo de viga i \in C ter� tempo de cura i.
%   limite_cura = tempo de cura maximo quando misturado for usado.
%   up_d = limite superior de demandas
function generate_inst(arq, C, M, up_d, misturado, V)
	clc;
	%M = 7;
	if misturado == 1
		%limite_cura = input('Digite um limite de tempo de cura:');
		limite_cura = 3;
	end
	
	%Dados usados para criar conjunto de formas
	caps = [11.95 5.95];
	%imprimir a amostra ordenada
	c = sort(datasample(caps, M));

	t = 0;
	soma_barras = 0;
	
	maior = 0;
	for i = 1:C
		%usar tempo de cura randomico de no maximo limite_cura se for
		%misturado
		if misturado == 1
			cura = ceil(rand(1,1)*limite_cura);
		else
			cura = i;
		end
		
		n_barras = ceil(rand(1,1)*3);
        if n_barras > maior
            maior = n_barras;
        end
		
		%dados usado para criar conjunto de tamanhos
		tams = [1.12 1.45 2.35 2.5 2.65 2.95 3.3];
        
        pegue_um_n = [6];
        peguei_um_n = datasample(pegue_um_n,1);
        
		tam = sort(datasample(tams,peguei_um_n,'Replace',false));

		N = length(tam);

		d = randi([ceil(up_d/3) up_d], 1, N);
		
		%estrutura usada para guardar dados do tipo gerado para a impress�o
		s(i) = struct('cura', {cura},'n_tam', {N}, 'tamanhos', {tam}, 'barras', {n_barras},...
			'demandas', {d}); %#ok<AGROW>
		%calculo do T
		t = cura * sum(tam.*d) + t;
		soma_barras = soma_barras + sum(tam.*d)*n_barras;
	end
	
	t = ceil(t/ sum(c));
	%aumento de 50% no T
	t = ceil(1.5*t);

	%gerar barras
	
	%Barras_tam = [12];
	Sobras_tam = [2 5 6 8];
	W = 1;
	Barras =  12;%datasample(Barras_tam, W, 'Replace',false);
	Sobras =  datasample(Sobras_tam, V, 'Replace',false);
    Sobras = sort(Sobras);
    
    UB_estoque = 2*t*M*maior;
    
    e = randi([ceil(UB_estoque/5) UB_estoque], 1, W+V);
    e(1) = UB_estoque;

    
    mkdir Instancia_novas
    arquivo =strcat('Instancias\',arq);
	%impress�o dos dados
	inst2 = fopen(arquivo, 'w');
	fprintf(inst2, '%d %d %d\n \n', C, M, t);
	
	fprintf(inst2, '%d %d \n \n', W, V);
	
	for i = 1:M
		fprintf(inst2, '%5.2f ', c(i));  
	end
	
	fprintf(inst2, '\n\n');
	
	for i = 1:C
		fprintf(inst2, '%d \n%d \n', s(i).cura, s(i).n_tam, s(i).barras);
		for j = 1:s(i).n_tam
			fprintf(inst2, '%5.2f ', s(i).tamanhos(j)); 
		end
		fprintf(inst2, '\n');
		for j = 1:s(i).n_tam
			fprintf(inst2, '%d ', s(i).demandas(j));  
		end
		fprintf(inst2, '\n\n');
	end
	
	for i = 1:W
		 fprintf(inst2, '%d ', Barras(i));
	end
	for i = 1:V
		 fprintf(inst2, '%d ', Sobras(i));
	end
	fprintf(inst2, '\n');
	for i = 1: W +V
		fprintf(inst2, '%d ', e(i));
	end
	fprintf(inst2, '\n%5.2f', 0.3);
	
	fclose(inst2);
end