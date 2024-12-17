clear all
clc

% Parâmetros do algoritmo memético
populacao_tamanho=6; % Tamanho da população
num_geracoes=100; % Número de gerações
prob_mutacao=0.5; % Probabilidade de mutação
margem_erro=0.11; % Margem de erro na avaliação da solução
distancia_minima=0.77; % Distância mínima entre soluções

% Inicialização
populacao=-120+(0-(-120))*rand(populacao_tamanho,1); % População inicial aleatória
melhores_solucoes=zeros(num_geracoes,1); % Vetor para armazenar os melhores indivíduos de cada geração
melhor_fitness=Inf;

% Função para verificar se uma solução viola as restrições
function viola=viola_restricoes(sol,solucoes,margem_erro,distancia_minima)
    viola=false;
    % Verifica se a solução viola a margem de erro
    if abs(sol)<margem_erro
        viola=true;
        return;
    end
    % Verifica a distância mínima em relação a outras soluções
    if any(abs(solucoes-sol)<distancia_minima)
        viola=true;
        return;
    end
end

% Função para realizar a busca local
function novo_x = busca_local(x, margem_erro)
    novo_x = x;
    f_x = f_vet(x);
    step_size = 0.1; % Tamanho do passo para a busca local
    max_iter = 100; % Número máximo de iterações
    for iter = 1:max_iter
        % Tentativa de mover para a esquerda
        x_left = x - step_size;
        f_left = f_vet(x_left);
        if abs(f_left) < abs(f_x) && abs(f_left) < margem_erro
            novo_x = x_left;
            f_x = f_left;
        else
            % Tentativa de mover para a direita
            x_right = x + step_size;
            f_right = f_vet(x_right);
            if abs(f_right) < abs(f_x) && abs(f_right) < margem_erro
                novo_x = x_right;
                f_x = f_right;
            else
                % Se nenhum movimento melhorar a solução, saia do loop
                break;
            end
        end
        x = novo_x;
    end
end

% Loop principal
for geracao=1:num_geracoes
    % Avaliação da população
    fitness=abs(f_vet(populacao));

    % Verifica se há uma solução melhor nesta geração
    [min_fitness,indice]=min(fitness);
    if min_fitness<melhor_fitness
        melhor_solucao=populacao(indice);
        melhor_fitness=min_fitness;
    end

    % Armazena o melhor indivíduo desta geração
    if geracao==1
        melhores_solucoes(geracao)=melhor_solucao;
    else
        % Verifica se o melhor indivíduo desta geração atende às restrições
        if ~viola_restricoes(melhor_solucao,melhores_solucoes(1:geracao-1),margem_erro,distancia_minima)
            melhores_solucoes(geracao)=melhor_solucao;
        end
    end

    % Seleção de pais (roleta simples)
    prob_selecao=1./fitness;
    prob_selecao=prob_selecao/sum(prob_selecao);
    pais_indices=randsample(1:populacao_tamanho,populacao_tamanho,true,prob_selecao);
    pais=populacao(pais_indices);

    % Reprodução (crossover)
    filhos=zeros(populacao_tamanho,1);
    for i=1:2:populacao_tamanho
        n=randi(2);
        pai1=pais(i);
        pai2=pais(i+1);
        alpha=rand();
        filhos(i)=alpha*pai1+(1-alpha)*((-1)^n)*pai2;
        filhos(i+1)=(1-alpha)*((-1)^n)*pai1+alpha*pai2;
    end

    % Mutação
    for i=1:populacao_tamanho
        if rand()<prob_mutacao
            filhos(i)=-120+(0-(-120))*rand();
        end
    end

    % Busca local
    for i=1:populacao_tamanho
        filhos(i) = busca_local(filhos(i), margem_erro);
    end

    % Atualiza a população
    populacao=filhos;
end

% Mostra o resultado
disp('Melhores valores de x encontrados em cada geração:');
disp(melhores_solucoes);

