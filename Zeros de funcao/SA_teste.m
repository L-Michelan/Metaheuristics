clear all
clc
% Parâmetros do algoritmo
T = 100; % Temperatura inicial
Tmin = 1; % Temperatura mínima
alpha = 0.99; % Taxa de resfriamento
iteracoes = 10; % Número de iterações

% Definição da função alvo


% Inicialização
x_best = []; % Vetor solução

% Laço principal
while T > Tmin
    for i = 1:iteracoes
        % Gera uma solução vizinha
        x_atual = -120 + (0-(-120)) * rand(); % Valor inicial aleatório dentro do domínio

        % Verifica se a solução vizinha está dentro do domínio
        if x_atual < -120 || x_atual > 0
            continue;
        end

        % Calcula o valor da função para a solução atual
        f_atual = f(x_atual);

        % Aceita o vizinho se ele for melhor ou com uma probabilidade
        % decrescente se ele for pior
        for j = 1:iteracoes
            % Gera uma solução vizinha
            x_vizinho = -120 + (0-(-120)) * rand(); % Valor aleatório dentro do domínio

            % Calcula o valor da função para a solução vizinha
            f_vizinho = f(x_vizinho);

            % Aceita o vizinho se ele for melhor ou com uma probabilidade
            % decrescente se ele for pior
            if f_vizinho < f_atual || rand() < exp((f_atual - f_vizinho)/T)
                x_atual = x_vizinho;
                f_atual = f_vizinho;
            end

            % Verifica se a solução vizinha cruzou o eixo zero com um erro de 0.11
            if abs(f_vizinho) < 0.77
                % Verifica se o valor encontrado já está presente no vetor solução
                if ~any(abs(x_best - x_vizinho) < 0.11)
                    % Adiciona o valor encontrado ao vetor solução
                    x_best = [x_best, x_vizinho];
                end
            end
        end
    end

    % Atualiza a temperatura
    T = alpha * T;
end
sol=unique(x_best);
% Mostra os resultados
disp('Valores de x que cruzam o eixo zero com erro de 0.11:');
disp(sol);
