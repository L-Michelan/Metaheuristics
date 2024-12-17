clear all
clc

% Parâmetros do algoritmo genético
populacao_tamanho = [4:2:20]; % Tamanho da população
num_geracoes = [50:10:100]; % Número de gerações
prob_mutacao = 0.5 * rand(1,5); % Probabilidade de mutação
margem_erro = [0.01:0.05:1]; % Margem de erro na avaliação da solução
distancia_minima = [0.1:0.05:0.5]; % Distância mínima entre soluções

% Função alvo
function y = f_vet(x)
    y = arrayfun(@(x) 2*exp(x) - x*sin(x+3) - 3, x);
end

% Função para verificar se uma solução viola as restrições
function viola = viola_restricoes(sol, solucoes, margem_erro, distancia_minima, L)
    viola = false;
    % Verifica a distância mínima em relação a outras soluções
    if any(abs(solucoes - sol) < distancia_minima(L))
        viola = true;
        return;
    end
end

melhor_contagem = 0;
melhores_parametros = [];

for I = 1:length(populacao_tamanho)
    for J = 1:length(num_geracoes)
        for K = 1:length(prob_mutacao)
            for L = 1:length(distancia_minima)
              for M = 1:length(margem_erro)
                % Inicialização
                populacao = -120 + (0 - (-120)) * rand(populacao_tamanho(I),1); % População inicial aleatória
                populacao = max(-120, min(0, populacao)); % Garantir que os valores da população inicial estejam dentro do intervalo [-120, 0]
                melhores_solucoes = zeros(num_geracoes(J),1); % Vetor para armazenar os melhores indivíduos de cada geração
                melhor_fitness = Inf;
                contador = 0;

                %guardando os valores de parametros
                pop=populacao_tamanho(I);
                geracoes=num_geracoes(J);
                mutacao=prob_mutacao(K);
                distancia=distancia_minima(L);
                margem=margem_erro(M);

                tic;
                % Loop principal
                for geracao = 1:num_geracoes(J)
                    % Avaliação da população
                    fitness = f_vet(populacao);

                    % Verifica se há uma solução melhor nesta geração
                    [min_fitness, indice] = min(abs(fitness));
                    if abs(min_fitness) < abs(melhor_fitness)
                        melhor_solucao = populacao(indice);
                        melhor_fitness = min_fitness;
                    end

                    % Armazena o melhor indivíduo desta geração
                    if geracao == 1
                        if abs(fitness(indice)) < margem_erro
                            melhores_solucoes(geracao) = melhor_solucao;
                            contador = contador + 1;
                        end
                    else
                        % Verifica se o melhor indivíduo desta geração atende às restrições
                        if ~viola_restricoes(melhor_solucao, melhores_solucoes(1:geracao-1), margem_erro, distancia_minima, L) && abs(fitness(indice)) < margem_erro
                            melhores_solucoes(geracao) = melhor_solucao;
                            contador = contador + 1;
                        end
                    end

                    % Seleção de pais (roleta simples)
                    prob_selecao = 1./fitness;
                    prob_selecao = prob_selecao / sum(prob_selecao);
                    pais_indices = randsample(1:populacao_tamanho(I), populacao_tamanho(I), true, prob_selecao);
                    pais = populacao(pais_indices);

                    % Reprodução (crossover)
                    filhos = zeros(populacao_tamanho(I),1);
                    for i = 1:2:populacao_tamanho(I)
                        n = randi(2);
                        pai1 = pais(i);
                        pai2 = pais(i+1);

                        % Inicializa os filhos fora do intervalo para garantir a entrada no loop
                        filho1 = -121;
                        filho2 = -121;

                        % Repete o cruzamento até que os filhos estejam dentro do intervalo
                        while filho1 < -120 || filho1 > 0 || filho2 < -120 || filho2 > 0
                            alpha = rand();
                            filho1 = alpha * pai1 + (1 - alpha) * ((-1)^n) * pai2;
                            filho2 = (1 - alpha) * ((-1)^n) * pai1 + alpha * pai2;

                            % Garante que os valores permaneçam dentro do intervalo [-120, 0]
                            filho1 = max(-120, min(0, filho1));
                            filho2 = max(-120, min(0, filho2));
                        end

                        filhos(i) = filho1;
                        filhos(i+1) = filho2;
                    end

                    % Mutação
                    for i = 1:populacao_tamanho(I)
                        if rand() < prob_mutacao(K)
                            % Mutação
                            % Adicionar uma mudança dentro do intervalo [-10, 10]
                            mudanca = -10 + (10 - (-10)) * rand();
                            filhos(i) = populacao(i) + mudanca;

                            % Garantir que os valores permaneçam dentro do intervalo [-120, 0]
                            filhos(i) = max(-120, min(0, filhos(i)));
                        end
                    end

                    % Atualiza a população
                    populacao = filhos;
                end
                tempo=toc;
                zeros_por_tempo=contador/tempo;

                %output de cada iteracao
                disp(['Tamanho da populacao ', num2str(pop), ', Número de geracoes: ', num2str(geracoes), ', Mutação: ', num2str(mutacao), ', Distância mínima: ', num2str(distancia), ', Erro: ', num2str(margem)]);

                disp(['Número de soluções ótimas encontradas: ', num2str(contador)]);
                % Verifica se o número de melhores soluções encontradas nesta iteração é maior que o máximo encontrado até agora
                if zeros_por_tempo > melhor_contagem
                    melhores_parametros = [populacao_tamanho(I), num_geracoes(J), prob_mutacao(K), distancia_minima(L), margem_erro(M)];
                    count = contador;
                    time=tempo;
                    melhor_contagem=zeros_por_tempo;
                end
                end
            end
        end
    end
end

% Exibe apenas os valores dos parâmetros do loop principal
disp(['melhor resultado obtido com:'])
melhores_parametros
