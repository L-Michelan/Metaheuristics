clear all
clc

%-----parametros e variaveis-----%
max_iter = 101; % numero maximo de iteracoes
Y = zeros(max_iter+1);% variavel para o grafico e resultado
X = 1:1:max_iter+1;

J = 16; % total de regioes da cidade a serem visitados
K = 2; % numero de drones
m = 8; % numero de locais potenciais para a instalacao das bases
dj = 400; % distancias para sobrevoo das regioes (metros)
Dk = [1860, 2580];  % capacidade de voo associadas a cada drone (em segundos)
tempo_recarga = [2760, 4200]; % tempo de recarga dos drones (em segundos)
velocidade = [4.72, 8.95]; % velocidade dos drones (em m/s)
dij = [10, 150, 309, 172, 285, 305, 435, 433; % distancias entre bases e regioes (metros)
       95, 10, 175, 172, 195, 320, 407, 350;
       230, 10, 24, 233, 168, 385, 407, 305;
       382, 170, 24, 357, 198, 492, 443, 305;
       10, 173, 309, 27, 230, 155, 300, 347;
       93, 83, 175, 27, 95, 282, 259, 235;
       230, 83, 24, 160, 20, 320, 259, 157;
       382, 190, 24, 312, 100, 417, 313, 157;
       165, 280, 350, 27, 230, 10, 197, 311;
       190, 230, 242, 27, 95, 95, 125, 173;
       285, 230, 170, 160, 10, 230, 125, 10;
       421, 290, 170, 312, 100, 385, 215, 10;
       303, 400, 435, 153, 275, 10, 150, 311;
       320, 368, 352, 153, 180, 95, 10, 175;
       382, 368, 305, 220, 153, 230, 10, 10;
       491, 407, 305, 350, 185, 385, 175, 10];
dij = dij';
djj = [ 0, 24, 148, 272, 24, 65, 190, 345, 148, 205, 272, 385, 272, 340, 385, 480; % distancia entre locais (metros)
        24, 0, 24, 148, 65, 24, 65, 210, 242, 148, 205, 272, 340, 272, 340, 385;
        148, 24, 0, 24, 190, 65, 24, 65, 272, 205, 148, 205, 385, 340, 272, 340;
        272, 148, 24, 0, 345, 190, 65, 24, 385, 272, 205, 148, 480, 385, 340, 272;
        24, 65, 190, 345, 0, 24, 148, 272, 24, 65, 190, 345, 148, 205, 272, 385;
        65, 24, 65, 210, 24, 0, 24, 148, 65, 24, 65, 190, 205, 148, 205, 272;
        190, 65, 24, 65, 148, 24, 0, 24, 190, 65, 24, 65, 272, 205, 148, 205;
        345, 210, 65, 24, 272, 148, 24, 0, 345, 190, 65, 24, 385, 272, 205, 148;
        148, 242, 272, 385, 24, 65, 190, 345, 0, 24, 148, 272, 24, 65, 190, 345;
        205, 148, 205, 272, 65, 24, 65, 190, 24, 0, 24, 148, 65, 24, 65, 190;
        272, 205, 148, 205, 190, 65, 24, 65, 148, 24, 0, 24, 190, 65, 24, 65;
        385, 272, 205, 148, 345, 190, 65, 24, 272, 148, 24, 0, 345, 190, 65, 24;
        272, 340, 385, 480, 148, 205, 272, 385, 24, 65, 190, 345, 0, 24, 148, 272;
        340, 272, 340, 385, 205, 148, 205, 272, 65, 24, 65, 190, 24, 0, 24, 148;
        385, 340, 272, 340, 272, 205, 148, 205, 65, 24, 65, 190, 65, 24, 0, 24;
        480, 385, 340, 272, 385, 272, 205, 148, 345, 190, 65, 24, 272, 148, 24, 0];
djj = djj;
%--------------------------------%


%-----inicializacao aleatoria-----%
function Xijk = inicializarAleatorio(m, J, K)
  Xijk = zeros(m, J, K);
  regioes = randperm(J); % gera um vetor de J numeros inteiros aleatorios únicos de 1 a J
  k = 1;
  for M = 1:J
    base(M) = randi(m); % selecionar aleatoriamente uma base para cada drone k
  endfor
  for j = 1:J
    Xijk(base(j), regioes(j), k) = 1;
    k = mod(k, K) + 1; % alterna entre os drones
  endfor
end
%---------------------------------%


%-----funcao objetivo-----%
function [f, f_individual] = calcularObjetivo(Xijk, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade)
  f_individual = zeros(1,K);
  f = 0;
  for k = 1:K
    [base, regiao] = find(Xijk(:,:,k) ~= 0); % encontra a base e a regiao associadas ao drone k
    retorno = false;
    for j = 1:length(regiao)
      if j == 1
        f_individual(k) = f_individual(k) + ((dij(base(1),regiao(j)) + dj) / velocidade(k)); % calculo do tempo gasto para o drone k ir da base i ate o regiao j
        Dk_bateria(k) = Dk_bateria(k) - (((dij(base(1),regiao(j)) + dj) / velocidade(k)) / Dk(k)); % desconta o tempo gasto da bateria restante do drone
      elseif Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dij(base(j),regiao(j)) + dj) / velocidade(k)) / Dk(k)) >= 0 % garante restricao 3, se resta tempo de bateria no drone para ir ate a proxima regiao, sobrevoa-la e voltar para a base (necessario somar a distancia ate a base para caso acabe a bateria do drone este nao fique ilhado)
        if retorno == true % se teve que voltar para a base
          f_individual(k) = f_individual(k) + ((dij(base(j - 1),regiao(j)) + dj) / velocidade(k));
          Dk_bateria(k) = Dk_bateria(k) - (((dij(base(j -1 ),regiao(j)) + dj) / velocidade(k)) / Dk(k));
          retorno = false;
        else
        f_individual(k) = f_individual(k) + ((djj(regiao(j - 1),regiao(j)) + dj) / velocidade(k)); % calculo do tempo gasto para o drone k ir da regiao j ate a proxima regiao j'
        Dk_bateria(k) = Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dj) / velocidade(k)) / Dk(k)); % desconta o tempo gasto da bateria restante do drone
        endif
      elseif j == length(regiao) && Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dij(base(j),regiao(j)) + dj) / velocidade(k)) / Dk(k)) >= 0 % se visitou todas as regioes j
        f_individual(k) = f_individual(k) + ((dij(base(j),regiao(j)) + dj + djj(regiao(j - 1),regiao(j))) / velocidade(k));
      else
        f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)) + (tempo_recarga(k) * (1 - Dk_bateria(k))); % calculo do tempo gasto para o drone voltar da regiao atual j ate a base i somado com o tempo de recarga do drone
        Dk_bateria(k) = Dk(k) / Dk(k); % recarrega o drone
        retorno = true;
      endif
    endfor
  endfor
  f = sum(f_individual(:));
end
%-------------------------%


%-----verificar restricoes-----%
function violacao = verificarRestricoes(X_ijk, J)
  violacao = false;
  for j = 1:J
    total_visitas = sum(sum(X_ijk(:,j,:))); % cada regiao j deve ser visitada por exatamente um drone uma unica vez
    if total_visitas ~= 1
      violacao = true;
      fprintf('Violação na região %d: total de visitas é %d\n', j, total_visitas);
      return
    endif
  endfor
  if sum(X_ijk(:)) ~= J % todas as regioes devem ser visitados
    violacao = true;
    return;
  endif
end
%------------------------------%


%-----busca local (troca os bairros entre os dois drones)-----%
function [Xijk] = buscaLocal(Xijk, K)
  % seleciona dois drones aleatórios
  drone1 = randi(K);
  drone2 = randi(K);
  while drone2 == drone1
    drone2 = randi(K); % garante que os drones selecionados sejam diferentes
  endwhile
  % obtem as regioes visitados pelos drones

  [linhas_d1, colunas_d1] = find(Xijk(:,:,drone1) ~= 0);
  [linhas_d2, colunas_d2] = find(Xijk(:,:,drone2) ~= 0);
  if length(colunas_d1) >= 2 && length(colunas_d2) >= 2 % seleciona aleatoriamente duas regioes para trocar entre os drones
    indices_d1 = randperm(length(colunas_d1), 2);
    indices_d2 = randperm(length(colunas_d2), 2);

    % localizacoes das regioes selecionados
    regioes_d1_troca = colunas_d1(indices_d1);
    regioes_d2_troca = colunas_d2(indices_d2);

    % troca as regioes entre os drones
    Xijk(linhas_d1(indices_d1(1)), regioes_d1_troca(1), drone1) = 0;
    Xijk(linhas_d1(indices_d1(2)), regioes_d1_troca(2), drone1) = 0;
    Xijk(linhas_d2(indices_d2(1)), regioes_d2_troca(1), drone2) = 0;
    Xijk(linhas_d2(indices_d2(2)), regioes_d2_troca(2), drone2) = 0;

    Xijk(linhas_d1(indices_d1(1)), regioes_d1_troca(1), drone2) = 1;
    Xijk(linhas_d1(indices_d1(2)), regioes_d1_troca(2), drone2) = 1;
    Xijk(linhas_d2(indices_d2(1)), regioes_d2_troca(1), drone1) = 1;
    Xijk(linhas_d2(indices_d2(2)), regioes_d2_troca(2), drone1) = 1;

endif
end
%-------------------------------------------------------%


%-----SA-----%
function [X_best, f_best, f_individual_best, Y] = SA(dj, Dk, m, J, K, max_iter, dij, djj, tempo_recarga, velocidade)
  T0 = 100; % temperatura inicial
  Tf = 0.1; % temperatura final
  alpha = 0.95; % taxa de resfriamento
  Dk_bateria = Dk ./ Dk; % inicia a bateria (iniciada como 1 para representar 100%)
  Xijk = inicializarAleatorio(m, J, K);
  while verificarRestricoes(Xijk, J)
    Xijk = inicializarAleatorio(m, J, K);
  endwhile
  X_best = Xijk;
  [f_best, f_individual_best] = calcularObjetivo(Xijk, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade)
  f_curr = 0;
  delta_f = f_curr - f_best;
  T = T0;
  iter = 1;
  busca = 0;
  mudanca = 0;
  aceito = false;
  Y(1) = f_best;

  while T > Tf && iter <= max_iter
    iter = iter + 1;
    if aceito == true
      X_ijk = X_aceito;
    else
      X_ijk = X_best;
    endif
    aceito = false;
    % aplica busca local (explotacao)
    X_ijk = buscaLocal(X_ijk, K);

     [f_curr, f_individual_curr] = calcularObjetivo(X_ijk, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade);
     delta_f = f_curr - f_best;

     if delta_f < 0 && ~verificarRestricoes(X_ijk, J)
       iter
        X_best = X_ijk;
        f_best = f_curr
     elseif exp(-delta_f / T) > rand() && ~verificarRestricoes(X_ijk, J)
        X_aceito = X_ijk;
        f_aceito = f_curr;
        aceito = true;
     endif
    T = alpha * T;
    Y(iter) = f_best;
  endwhile
end
%------------%

% chamar o algoritmo SA
tic
[X_best, f_best, f_individual_best, Y] = SA(dj, Dk, m, J, K, max_iter, dij, djj, tempo_recarga, velocidade);
tempo_CPU = toc
% exibir resultado
disp('Melhor solução encontrada:');
disp(X_best);
disp(['Tempo total: ', num2str(f_best)]);
disp(['Tempo de cada drone: ', num2str(f_individual_best)]);

plot(X,Y);
title('Curva de convergência SA');
xlabel('Iteração');
ylabel('Valor obtido');