function [best_solution, best_cost, initial_solution, best_eq_values] = vns_sistema_nao_linear()

    % Define os parâmetros do algoritmo
    max_iter = 1000; % Número máximo de iterações
    neighborhood_sizes = [1, 2, 3]; % Tamanhos das vizinhanças
    neighborhood_iter = 50; % Número de iterações em cada vizinhança

    % Define os limites para as variáveis
    lower_bound = [0 0 0 0 0.1 0.1 0 0 0 0 0 0 0]; % Limites inferiores ajustados para x5 e x6 >= 0.1 e x7 a x13 >= 0
    upper_bound = [100 100 100 100 10 10 10 10 10 10 10 10 10]; % Limites superiores ajustados para x1 a x4 = 100, x5 e x6 >= 0.1 e x7 a x13 >= 0

    % Gera uma solução inicial aleatória dentro dos limites
    initial_solution = generate_initial_solution(lower_bound, upper_bound);
    best_solution = initial_solution;

    % Calcula o custo da solução inicial
    [best_cost, best_eq_values] = cost_function(best_solution);

    % Inicia o loop principal do algoritmo
    for iter = 1:max_iter
        % Para cada tamanho de vizinhança
        for neighborhood_size = neighborhood_sizes
            % Gera uma solução vizinha aleatória dentro da vizinhança
            neighbor_solution = generate_neighbor(best_solution, neighborhood_size, lower_bound, upper_bound);

            % Calcula o custo da solução vizinha
            [neighbor_cost, neighbor_eq_values] = cost_function(neighbor_solution);

            % Aceita a solução vizinha se ela for melhor
            if neighbor_cost < best_cost
                best_solution = neighbor_solution;
                best_cost = neighbor_cost;
                best_eq_values = neighbor_eq_values;
                break; % Interrompe o loop sobre os tamanhos de vizinhança após encontrar uma melhoria
            end
        end
    end

    % Retorne os valores iniciais e os melhores valores encontrados
    disp('Melhor solução encontrada:');
    disp(best_solution);
    disp('Custo da melhor solução:');
    disp(best_cost);
    disp('Solução inicial:');
    disp(initial_solution);
    disp('Valores das equações para a melhor solução encontrada:');
    disp(best_eq_values);
end

function initial_solution = generate_initial_solution(lower_bound, upper_bound)
    % Gera uma solução inicial aleatória dentro dos limites
    initial_solution = zeros(size(lower_bound)); % Use o tamanho de lower_bound para garantir consistência
    initial_solution(1:4) = randi([1, 100], 1, 4); % Atribui valores inteiros positivos para x1 a x4 com limite de 100
    initial_solution(5:end) = lower_bound(5:end) + randn(size(lower_bound(5:end))).*0.1 .* (upper_bound(5:end) - lower_bound(5:end)); % Perturbação gaussiana para x5 a x13
    initial_solution(13) = 0; % Definindo x13
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

function neighbor = generate_neighbor(current_solution, neighborhood_size, lower_bound, upper_bound)
    % Gera uma solução vizinha aleatória dentro da vizinhança
    neighbor = current_solution;
    for i = 1:neighborhood_size
        % Perturbação gaussiana para x5 a x13
        neighbor(5:end) = neighbor(5:end) + randn(size(neighbor(5:end))).*0.1 .* (upper_bound(5:end) - lower_bound(5:end));
        % Perturbação discreta para x1 a x4
        neighbor(1:4) = max(1, round(neighbor(1:4) + randn(size(neighbor(1:4))))); % Garante que x1 a x4 sejam inteiros positivos
        % Garante que a vizinhança está dentro dos limites
        neighbor = max(neighbor, lower_bound);
        neighbor = min(neighbor, upper_bound);
    end
end

% Executa o algoritmo Variable Neighborhood Search
[best_solution, best_cost, initial_solution, best_eq_values] = variable_neighborhood_search();

