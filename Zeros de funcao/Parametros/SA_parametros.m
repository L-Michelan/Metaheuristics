clear all
clc

% Definindo a função alvo
function y = f(x)
    y = 2*exp(x) - x*sin(x+3) - 3;
end

% Parâmetros do algoritmo
T = [100:10:200]; % Temperatura inicial
Tmin = 1; % Temperatura mínima
alpha = (0.95 - 0.5) * rand(5, 1) + 0.5; % Taxa de resfriamento
iteracoes = [10:10:100]; % Número de iterações
distancia_minima = [0.1:0.05:0.5];
margem_erro = [0.01:0.05:1];
desvios_padrao = [0.1:0.1:1]; % Desvios padrão para teste

% Inicialização
x_best = {}; % Vetor solução
melhor_contagem = 0;
melhores_parametros = {}; % Inicializa a variável como um vetor vazio
zeros_por_iteracao=0;
for I = 1:length(iteracoes)
    for J = 1:length(T)
        for K = 1:length(margem_erro)
            for O = 1:length(distancia_minima)
                for L = 1:length(alpha)
                    for M = 1:length(desvios_padrao)
                        % Parâmetros atuais
                        iter = iteracoes(I);
                        temp = T(J);
                        erro = margem_erro(K);
                        dist_min = distancia_minima(O);
                        alpha_atual = alpha(L);
                        desvio_padrao = desvios_padrao(M);

                        contador = 0; % Inicializa o contador para cada conjunto de parâmetros
                        x_best = {}; % Vetor solução

                        % Laço principal do Simulated Annealing
                        while temp > Tmin
                            % Inicializa x_atual com um valor aleatório dentro do intervalo [-120, 0]
                            x_atual = (-120 - 0) * rand() + 0;

                            % Gera uma solução vizinha única
                            x_vizinho = x_atual + randn() * desvio_padrao;
                            while x_vizinho < -120 || x_vizinho > 0 || any(abs(cellfun(@(x) norm(x - x_vizinho), x_best)) < dist_min)
                                x_vizinho = x_atual + randn() * desvio_padrao;
                            end

                            % Calcula o valor da função para a solução vizinha
                            f_vizinho = f(x_vizinho);

                            % Calcula a variação de energia
                            deltaE = f_vizinho - f(x_atual);

                            % Aceita a solução vizinha se ela for melhor ou com probabilidade e^(-deltaE/temp)
                            if deltaE <= 0 || rand() < exp(-deltaE / temp)
                                x_atual = x_vizinho;
                                % Verifica se a solução atende aos critérios de erro e distância mínima
                                if abs(f_vizinho) <= erro
                                    x_best{end+1} = x_vizinho;
                                    contador = contador + 1;
                                    zeros_por_iteracao=contador/iter;
                                end
                            end

                            % Atualiza a temperatura
                            temp = alpha_atual * temp;
                        end % while temp > Tmin

                        % Output com os parâmetros utilizados na iteração atual do while loop
                        disp(['Iterações: ', num2str(iter), ', Temperatura: ', num2str(T(J)), ', Erro: ', num2str(erro), ', Distância mínima: ', num2str(dist_min), ', Alpha: ', num2str(alpha_atual), ', Desvio Padrão: ', num2str(desvio_padrao)]);

                        % Exibe o número de vezes que encontrou uma solução x_best e as soluções encontradas
                        disp(['Número de soluções ótimas encontradas: ', num2str(contador)]);
                        disp('Soluções encontradas de x_best:');
                        disp([x_best{:}]);

                        % Verifica se o contador atual é melhor do que o melhor contador encontrado até agora
                        if zeros_por_iteracao > melhor_contagem
                            melhores_parametros = [iter, temp, erro, dist_min, desvio_padrao];
                            melhor_contagem = zeros_por_iteracao;
                        end
                    end
                end
            end
        end
    end
end

% Verifica se algum conjunto de parâmetros atendeu aos critérios
if ~isempty(melhores_parametros)
    disp(['Número de iterações: ', num2str(melhores_parametros(1)), ', Temperatura inicial: ', num2str(melhores_parametros(2)), ', Margem de erro: ', num2str(melhores_parametros(3)), ', Distância mínima: ', num2str(melhores_parametros(4)),'Desvio Padrão: ', num2str(melhores_parametros(5)), ', Número de soluções ótimas encontradas: ', num2str(melhor_contagem)]);
else
    disp('Nenhum conjunto de parâmetros atendeu aos critérios.');
end

