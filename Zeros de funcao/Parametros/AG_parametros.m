% Função objetivo modificada para contar cruzamentos com o eixo x=0
function n_crossings = f_vet(x)
    y = 2*exp(x) - x*sin(x+3) - 3;
    n_crossings = numel(find(y(1:end-1).*y(2:end) <= 0));
end

% Função para avaliar a aptidão (fitness) de uma solução com base nos cruzamentos
function fitness = avaliar_aptidao(solucoes)
    fitness = f_vet(solucoes);
end

% Algoritmo Genético para ajustar os parâmetros
function parametros_ajustados = ajustar_parametros_AG(populacao_tamanho, num_geracoes, prob_mutacao, margem_erro, distancia_minima)
    % Inicialização
    populacao = -120 + (0 - (-120)) * rand(populacao_tamanho, 1); % População inicial de parâmetros aleatórios

    % Loop principal do algoritmo genético
    for geracao = 1:num_geracoes
        % Avaliação da aptidão (fitness) de cada indivíduo na população
        fitness = avaliar_aptidao(populacao);

        % Seleção dos melhores indivíduos (elitismo)
        [~, indice] = max(fitness);
        melhor_parametro = populacao(indice);

        % Atualização dos melhores parâmetros encontrados
        if geracao == 1 || f_vet(melhor_parametro) > f_vet(melhor_parametros_ajustados)
            melhor_parametros_ajustados = melhor_parametro;
        end

        % Crossover e mutação (opcional)
        % Aqui você pode implementar operadores genéticos como crossover e mutação

        % Exemplo de mutação: adicionar um ruído aleatório aos parâmetros
        populacao = populacao + randn(populacao_tamanho, 1) * 0.1;

        % Limitando os parâmetros dentro de limites razoáveis (opcional)
        % Aqui você pode aplicar limites aos parâmetros, se necessário

        % Exemplo de limitação dos parâmetros:
        populacao = max(populacao, -120); % Limitando a -120

        % Exemplo de limitação dos parâmetros:
        populacao = min(populacao, 0); % Limitando a 0
    end

    % Retorno dos parâmetros ajustados
    parametros_ajustados = struct('populacao_tamanho', populacao_tamanho, ...
                                  'num_geracoes', num_geracoes, ...
                                  'prob_mutacao', prob_mutacao, ...
                                  'margem_erro', margem_erro, ...
                                  'distancia_minima', distancia_minima, ...
                                  'num_cruzamentos', f_vet(melhor_parametros_ajustados));
end

% Parâmetros do algoritmo genético
populacao_tamanho = 6; % Tamanho da população
num_geracoes = 100; % Número de gerações
prob_mutacao = 0.5; % Probabilidade de mutação
margem_erro = 0.11; % Margem de erro na avaliação da solução
distancia_minima = 0.77; % Distância mínima entre soluções

% Ajuste dos parâmetros do algoritmo genético
parametros_ajustados = ajustar_parametros_AG(populacao_tamanho, num_geracoes, prob_mutacao, margem_erro, distancia_minima);

% Mostra o resultado
disp('Parâmetros ajustados:');
disp(parametros_ajustados);
