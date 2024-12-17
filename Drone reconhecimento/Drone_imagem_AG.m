clear all
clc

%-----parametros e variaveis-----%
geracoes = 100; % numero maximo de iteracoes
pop_size = 10; % tamanho da populacao

Y = zeros(1, geracoes+1);% variavel para o grafico e resultado
X = 1:1:geracoes+1;

num_pares = pop_size / 2; % define o número de pares de pais (5 neste caso)
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
  regioes = randperm(J); % gera um vetor de J numeros inteiros aleatorios unicos de 1 a J

  % determina quantas regioes cada drone deve cobrir para garantir que um drone nao fique muito sobrecarregado
  regioes_por_drone = J / K;

  for k = 1:K
    % seleciona as regioes atribuidas ao drone k
    regioes_drone = regioes((k-1)*regioes_por_drone + 1:k*regioes_por_drone);

    % seleciona aleatoriamente uma base para o drone k
    base = randi(m);

    % associa as regioes selecionadas ao drone k
    for j = regioes_drone
      Xijk(base, j, k) = 1;
    endfor
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
          Dk_bateria(k) = Dk_bateria(k) - (((dij(base(j - 1),regiao(j)) + dj) / velocidade(k)) / Dk(k));
          retorno = false
        else
        f_individual(k) = f_individual(k) + ((djj(regiao(j - 1),regiao(j)) + dj) / velocidade(k)); % calculo do tempo gasto para o drone k ir da regiao j ate a proxima regiao j'
        Dk_bateria(k) = Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dj) / velocidade(k)) / Dk(k)); % desconta o tempo gasto da bateria restante do drone
        endif
      elseif j == length(regiao) && Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dij(base(1),regiao(j)) + dj) / velocidade(k)) / Dk(k)) >= 0 % se visitou todas as regioes j
        f_individual(k) = f_individual(k) + ((dij(base(j),regiao(j)) + dj + djj(regiao(j - 1),regiao(j))) / velocidade(k));
      else
        f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)) + (tempo_recarga(k) * (1 - Dk_bateria(k))); % calculo do tempo gasto para o drone voltar da regiao atual j ate a base i somado com o tempo de recarga do drone
        Dk_bateria(k) = Dk(k) / Dk(k); % recarrega o drone
        retorno = true
      endif
    endfor
  endfor
  f = sum(f_individual(:));
end
%---------------------------------%


%-----verificar restricoes-----%
function violacao = verificarRestricoes(X_ijk, J)
  violacao = false;
  for j = 1:J
    total_visitas = sum(sum(X_ijk(:,j,:))); % cada regiao j deve ser visitada por exatamente um drone uma unica vez
    if total_visitas ~= 1
      violacao = true;
      fprintf('Violação na região %d: total de visitas é %d\n', j, total_visitas);
      return;
    endif
  endfor
  if sum(X_ijk(:)) ~= J % todas as regioes devem ser visitados
    violacao = true;
    return;
  endif

  % um drone devem visitar ao menos 6 regioes para nao sobrecarregar o outro drone
  if sum(sum(X_ijk(:,:,1))) < 6
    violacao = true;
    return;
  elseif sum(sum(X_ijk(:,:,2))) < 6
    violacao = true;
    return;
  endif
end

% checa se os cada drone eh responsavel por pelo menos 45% da carga de trabalho, ou seja, visita pelo menos 20% dos locais infectados (45% para ter uma margem para o código nao violar as restricoes sempre que tiver um resultado nao perfeito)
function violacao2 = restricao_sobrecarga(X_ijk, J)
  violacao2 = false;
  if sum(sum(X_ijk(:,:,1))) < floor(0.45 * J) || sum(sum(X_ijk(:,:,2))) < floor(0.45 * J)
    violacao2 = true;
    return;
  endif
end
%------------------------------%



%-----crossover-----%
function filho = crossover(pai1, pai2, m, J, K)
  ponto_corte = randi([1 J-1]); % ponto de corte aleatorio
  filho = zeros(m, J, K); % inicializa o filho
  for k = 1:K
    filho(:,1:ponto_corte,k) = pai1(:,1:ponto_corte,k);
    filho(:,ponto_corte+1:end,k) = pai2(:,ponto_corte+1:end,k);
  endfor
end
%--------------------%


%-----mutacao (seleciona aleatoriamente 2 novas bases)-----%
function individuo_mutado = mutacao(individuo, m, J, K, probabilidade_mutacao)
  individuo_mutado = individuo;
  if rand() < probabilidade_mutacao
    % implementacao da mutacao de troca aleatoria de bases
    idx = randi([1 m], 1, 2); % gera dois indices aleatorios para bases
    temp = individuo_mutado(idx(1), :, :); % armazena temporariamente a base escolhida
    individuo_mutado(idx(1), :, :) = individuo_mutado(idx(2), :, :); % troca as bases
    individuo_mutado(idx(2), :, :) = temp; % finaliza a troca
  endif
end
%-------------------%


%---selecao de pais para cruzamento---%
function pais = selecao(populacao, f, num_pares)
  % ordena a populacao com base na funcao objetivo
  [~, indices] = sort(f);

  % seleciona os melhores individuos da populacao
  selecionados = indices(1:num_pares*2); % escolhe 5 pares de pais

  % cria matriz para armazenar os pares de pais
  pais = zeros(size(populacao, 1), size(populacao, 2), size(populacao, 3), num_pares*2);

  for i = 1:num_pares*2
    pais(:,:,:,i) = populacao(:,:,:,selecionados(i));
  endfor
end
%-------------------------------------%


%-----algoritmo genetico-----%
tic
probabilidade_mutacao = 0.05;

% inicializa a populacao
populacao = zeros(m, J, K, pop_size);
f = zeros(pop_size, 1);

for i = 1:pop_size
  populacao(:,:,:,i) = inicializarAleatorio(m, J, K);
endfor

% loop de geracoes
for iter = 1:geracoes
  Dk_bateria = Dk ./ Dk;
  % avalia a populacao
  for i = 1:pop_size
    % vrifica restricoes
    if verificarRestricoes(populacao(:,:,:,i), J)
      f(i) = Inf; % penaliza solucoes que violam as restricoes
    elseif restricao_sobrecarga(populacao(:,:,:,i), J)
      [f(i), ~] = calcularObjetivo(populacao(:,:,:,i), dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
      f(i) += 600;
    else
      [f(i), ~] = calcularObjetivo(populacao(:,:,:,i), dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
    endif
  endfor
  if iter == 1
    inicial = 1
    min(f)
    best_fitness(iter) = min(f);
    Y(iter) = best_fitness(iter);
    atual = min(f);
  elseif min(f) < atual
    iter
    min(f)
    best_fitness(iter + 1) = min(f);
    atual = min(f);
  else
    best_fitness(iter + 1) = best_fitness(iter);
  endif

% atualizacao da chamada da funcao de selecao
pais = selecao(populacao, f, num_pares);

% crossover e mutacao
  for i = 1:2:pop_size
    pai1 = pais(:,:,:,i);
    pai2 = pais(:,:,:,i+1);
    sobrecarga = true;
    count = 0;
    while sobrecarga == true
      sobrecarga = false;
      count = count + 1;
      filho1 = crossover(pai1, pai2, m, J, K);
      filho2 = crossover(pai2, pai1, m, J, K);
      for r = 1:K
        if (sum(sum(filho1(:,:,r))) < floor(0.45 * J) && count < 10) || (sum(sum(filho2(:,:,r))) < floor(0.45 * J) && count < 10)
          sobrecarga = true;
        elseif count == 10
          sobrecarga = false;
        endif
      endfor
    endwhile

      % avalia novos individuos
    if verificarRestricoes(filho1, J)
      f1 = Inf; % penaliza solucoes que violam as restricoes
    elseif restricao_sobrecarga(filho1)
      [f1, ~] = calcularObjetivo(filho1, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
      f1 += 600;
    else
      [f1, ~] = calcularObjetivo(filho1, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
    endif

    if verificarRestricoes(filho2, J)
      f2 = Inf; % penaliza solucoes que violam as restricoes
    elseif restricao_sobrecarga(filho2)
      [f2, ~] = calcularObjetivo(filho2, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
      f2 += 600;
    else
      [f2, ~] = calcularObjetivo(filho2, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
    endif
    if f1 < atual && f1 < f2
      iter
      f1
      best_fitness(iter + 1) = f1;
      atual = f1;
    elseif f2 < atual && f2 < f1
      iter
      f2
      best_fitness(iter + 1) = f2;
      atual = f2;
    endif

  % alicar mutacao se chance for cumprida
    if probabilidade_mutacao > rand
      filho_mutado1 = mutacao(filho1, m, J, K, probabilidade_mutacao);
      filho_mutado2 = mutacao(filho2, m, J, K, probabilidade_mutacao);

      % avalia novos individuos
      if verificarRestricoes(filho_mutado1, J)
        f_mutado1 = Inf; % penaliza solucoes que violam as restricoes
      elseif restricao_sobrecarga(filho_mutado1)
        [f_mutado1, ~] = calcularObjetivo(filho_mutado1, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
        f_mutado1 += 600;
      else
        [f_mutado1, ~] = calcularObjetivo(filho_mutado1, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
      endif

      if verificarRestricoes(filho_mutado2, J)
        f_mutado2 = Inf; % penaliza solucoes que violam as restricoes
      elseif restricao_sobrecarga(filho_mutado2)
        [f_mutado2, ~] = calcularObjetivo(filho_mutado2, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
        f_mutado2 += 600;
      else
        [f_mutado2, ~] = calcularObjetivo(filho_mutado2, dj, K, dij, djj, Dk, Dk_bateria, tempo_recarga, velocidade);
      endif

      if f_mutado1 < atual && f_mutado1 < f_mutado2
        iter
        f_mutado1
        best_fitness(iter + 1) = f_mutado1;
        atual = f_mutado1;
      elseif f_mutado2 < atual && f_mutado2 < f_mutado1
        iter
        f_mutado2
        best_fitness(iter + 1) = f_mutado2;
        atual = f_mutado2;
      endif

    % seleciona os melhores filhos para a proxima populacao (entre individuo cruzado e individuo sob selecao cultural)
      filho1 = filho_mutado1;
      filho2 = filho_mutado2;
    endif
    populacao(:,:,:,i) = filho1;
    populacao(:,:,:,i+1) = filho2;
  endfor
  if iter == 1
    best_fitness(iter + 1) = atual;
  endif
  Y(iter + 1) = best_fitness(iter);
endfor

% exibir valores
Tempo = toc
disp(['Tempo total: ', num2str(min(best_fitness))]);

plot(X,Y);
title('Curva de convergência AG');
xlabel('Geração');
ylabel('Melhor fitness');
