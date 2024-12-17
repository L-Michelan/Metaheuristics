clear all
clc

function [best_solution, best_cost, best_eq_values, first_best_solution] = ga_sistema_nao_linear()

    % Define os parâmetros do algoritmo genético
    population_size = 50; % Tamanho da população
    max_generations = 100; % Número máximo de gerações
    mutation_rate = 0.1; % Taxa de mutação
    crossover_rate = 0.8; % Taxa de crossover

    % Define os limites para as variáveis
    lower_bound = [0 0 0 0 0.1 0.1 0 0 0 0 0 0 0]; % Limites inferiores ajustados para x5 e x6 >= 0.1 e x7 a x13 >= 0
    upper_bound = [100 100 100 100 10 10 10 10 10 10 10 10 10]; % Limites superiores ajustados para x1 a x4 = 100, x5 e x6 >= 0.1 e x7 a x13 >= 0

    % Inicializa a população aleatória dentro dos limites
    population = initialize_population(population_size, lower_bound, upper_bound);

    % Avalia a população inicial
    [costs, eq_values] = evaluate_population(population);

    % Encontra a primeira melhor solução
    [first_best_cost, first_best_index] = min(costs);
    first_best_solution = population(first_best_index, :);

    % Inicia o loop principal do algoritmo genético
    for generation = 1:max_generations
        % Seleciona os pais para o crossover
        parents = selection(population, costs);

        % Realiza o crossover para gerar filhos
        children = crossover(parents, crossover_rate);

        % Aplica mutação aos filhos
        mutated_children = mutation(children, mutation_rate, lower_bound, upper_bound);

        % Avalia a nova população de filhos mutados
        [mutated_costs, mutated_eq_values] = evaluate_population(mutated_children);

        % Substitui a população antiga pela nova população de filhos mutados, se forem melhores
        better_indices = find(mutated_costs < costs);
        population(better_indices, :) = mutated_children(better_indices, :);
        costs(better_indices) = mutated_costs(better_indices);
        eq_values(better_indices, :) = mutated_eq_values(better_indices, :);
    end

    % Encontra a melhor solução na população final
    [best_cost, best_index] = min(costs);
    best_solution = population(best_index, :);
    best_eq_values = eq_values(best_index, :);

    % Mostra os resultados
    disp('Melhor solução encontrada:');
    disp(best_solution);
    disp('Custo da melhor solução:');
    disp(best_cost);
    disp('Valores das equações para a melhor solução encontrada:');
    disp(best_eq_values);
    disp('Primeira melhor solução:');
    disp(first_best_solution);
    disp('Custo da primeira melhor solução:');
    disp(first_best_cost);
end

function population = initialize_population(population_size, lower_bound, upper_bound)
    % Inicializa a população aleatória dentro dos limites
    num_variables = numel(lower_bound);
    population = zeros(population_size, num_variables);
    for i = 1:population_size
        population(i, 1:4) = randi([1, 100], 1, 4); % Atribui valores inteiros positivos para x1 a x4 com limite de 100
        population(i, 5:end) = lower_bound(5:end) + rand(size(lower_bound(5:end))) .* (upper_bound(5:end) - lower_bound(5:end));
    end
end

function [costs, eq_values] = evaluate_population(population)
    % Avalia a população
    num_individuals = size(population, 1);
    costs = zeros(num_individuals, 1);
    eq_values = zeros(num_individuals, 5); % Existem 5 equações no total
    for i = 1:num_individuals
        [costs(i), eq_values(i, :)] = cost_function(population(i, :));
    end
end

function [parents] = selection(population, costs)
    % Seleciona os pais com base no torneio binário
    num_parents = size(population, 1);
    tournament_size = 2; % Tamanho do torneio
    parents = zeros(num_parents, size(population, 2));
    for i = 1:num_parents
        % Seleciona aleatoriamente dois indivíduos para o torneio
        tournament_indices = randperm(num_parents, tournament_size);
        [~, winner_index] = min(costs(tournament_indices));
        parents(i, :) = population(tournament_indices(winner_index), :);
    end
end

function [children] = crossover(parents, crossover_rate)
    % Realiza o crossover de ponto único
    num_parents = size(parents, 1);
    num_variables = size(parents, 2);
    num_children = num_parents;
    children = zeros(num_children, num_variables);
    for i = 1:2:num_parents
        if rand() < crossover_rate
            crossover_point = randi([1, num_variables - 1]);
            children(i, :) = [parents(i, 1:crossover_point), parents(i + 1, crossover_point + 1:end)];
            children(i + 1, :) = [parents(i + 1, 1:crossover_point), parents(i, crossover_point + 1:end)];
        else
            children(i, :) = parents(i, :);
            children(i + 1, :) = parents(i + 1, :);
        end
    end
end

function [mutated_children] = mutation(children, mutation_rate, lower_bound, upper_bound)
    % Aplica mutação por adição gaussiana
    num_children = size(children, 1);
    num_variables = size(children, 2);
    mutated_children = children;
    for i = 1:num_children
        for j = 1:num_variables
            if rand() < mutation_rate
                mutated_children(i, j) = children(i, j) + randn() * 0.1 * (upper_bound(j) - lower_bound(j));
            end
        end
    end
    % Garante que os valores mutados permaneçam dentro dos limites
    mutated_children = max(mutated_children, lower_bound);
    mutated_children = min(mutated_children, upper_bound);
end

function [cost, eq_values] = cost_function(x)
    % Sistema de equações
    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    x4 = x(4);
    x5 = x(5);
    x6 = x(6);
    x7 = x(7);
    x8 = x(8);
    x9 = x(9);
    x10 = x(10);
    x11 = x(11);
    x12 = x(12);
    x13 = x(13);

    A1 = -6*x1 - 5*x2 + 102;
    A2 = -2*x3 - x4 + 200;
    A3 = -x2 - 2*x3 - 3*x4 + 150;

    eq1 = (5 + 200*pi*sin(pi*x1)*cos(pi*x1))*x6 + x5 - 6*x7 - x10;
    eq2 = (5 + 200*pi*sin(pi*x2)*cos(pi*x2))*x6 + x5 - 5*x7 - x9 - x11;
    eq3 = (5 + 200*pi*sin(pi*x3)*cos(pi*x3))*x6 + x5 - 2*x8 - 2*x9 - x12;
    eq4 = (5 + 200*pi*sin(pi*x4)*cos(pi*x4))*x6 + x5 - x8 - 3*x9 - x13;
    eq5 = x7*A1 + x8*A2 + x9*A3 - x10*x1 - x11*x2 - x12*x3 - x13*x4;

    eq_values = [eq1, eq2, eq3, eq4, eq5]; % Salva os valores das equações
    cost = eq1^2 + eq2^2 + eq3^2 + eq4^2 + eq5^2;
end

% Executa o algoritmo genético
[best_solution, best_cost, best_eq_values, first_best_solution] = ga_sistema_nao_linear();

