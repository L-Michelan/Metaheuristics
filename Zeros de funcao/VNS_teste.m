clear all
clc

% Parâmetros do algoritmo
max_iteracoes=100; % Número máximo de iterações
alfa=0.11; % Tamanho do vizinho

% Definição da função alvo
f=@(x)2*exp(x)-x*sin(x+3)-3;

% Inicialização
x_best=[]; % Vetor solução

% Laço principal
for iter=1:max_iteracoes
    % Gera uma solução inicial aleatória dentro do domínio
    x_atual=-120+(0-(-120))*rand(); % Valor inicial aleatório dentro do domínio

    % Aplica o algoritmo VNS
    for i=1:max_iteracoes
        % Calcula o valor da função para a solução atual
        f_atual=f(x_atual);
        dominio=false;
            while dominio==false
            x_vizinho= xatual + alfa * (2 * rand() − 1); % Valor aleatório dentro do domínio
            if (x_vizinho>-120 && x_vizinho<0)
              dominio=true;
            else
              dominio=false;
            endif
          endwhile
        % Calcula o valor da função para a solução vizinha
        f_vizinho=f(x_vizinho);

        % Aceita o vizinho se ele for melhor
        if f_vizinho<f_atual
            x_atual=x_vizinho;
            f_atual=f_vizinho;
        endif

        % Verifica se a solução vizinha cruzou o eixo zero com um erro de 0.11
        if abs(f_vizinho)<0.77
            % Verifica se o valor encontrado já está presente no vetor solução
            if ~any(abs(x_best-x_vizinho)<0.11)
                % Adiciona o valor encontrado ao vetor solução
                x_best=[x_best,x_vizinho];
            endif
        endif
    endfor
endfor

% Remove duplicatas dos valores encontrados
x_best=unique(x_best);

% Mostra os resultados
disp('Valores de x que cruzam o eixo zero com erro de 0.11:');
disp(x_best);
