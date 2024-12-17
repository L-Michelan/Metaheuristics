clear all
clc

%-----parametros e variaveis-----%
pop_size = 10; % tamanho da populacao
geracoes = 100; % quantidade de geracoes
Y = zeros(1, geracoes+1);% variavel para o grafico e resultado
X = 1:1:geracoes+1;

J = 16; % total de regioes da cidade a serem visitados
K = 4; % numero de drones
m = 8; % numero de locais potenciais para a instalacao das bases
num_pares = pop_size / 2; % define o numero de pares de pais para a selecao (5 neste caso)

dj = 400; % area de sobrevoo das casas (metros quadrados)
Dk = [468, 468, 600, 600];  % capacidade de voo associadas a cada drone (em segundos)
Tk = [30, 30, 16, 16]; % capacidade do tanque de cada drone (em L)
Vk = [0.12, 0.12, 0.093, 0.093]; % taxa de vazao de liquido de cada drone (em L/s)
Vp = 28.9 * (10^-3); % taxa de aplicacao do produto (em L/metro quadrado)
tempo_recarga = [720, 720, 1800, 1800]; % tempo de recarga dos drones (em segundos)
velocidade = [7, 7, 12, 12]; % velocidade dos drones (em m/s)
ppc = Vp * dj; %quantidade de produto necessario por casa

dij = [10, 150, 309, 172, 285, 305, 435, 433; % distancias entre bases e locais (metros)
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

dij = dij'; % toda matriz que estao transposta serve apenas para que os indices da matriz referencie exatamente aos indices da imagem da cidade hipotetica
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

dcc = [0, 40, 80, 40, 56.57, 89.44, 80, 89.44, 113.14; % distancia entre casas de um bairro
       40, 0, 40, 56.57, 40, 56.57, 89.44, 80, 89.44;
       80, 40, 0, 89.44, 56.57, 40, 113.14, 89.44, 80;
       40, 56.57, 89.44, 0, 40, 80, 40, 56.57, 89.44;
       56.57, 40, 56.57, 40, 0, 40, 56.57, 40, 56.57;
       89.44, 56.57, 40, 80, 40, 0, 89.44, 56.57, 40;
       80, 89.44, 113.14, 40, 56.57, 89.44, 0, 40, 80;
       89.44, 80, 89.44, 56.57, 40, 56.57, 40, 0, 40;
       113.14, 89.44, 80, 89.44, 56.57, 40, 80, 40, 0];

% Define as matrizes de bairros e casas infectadas
bairros = [1 0 0 1; 1 1 0 1; 0 1 0 1; 0 0 1 1];
indice = find(bairros'~=0);

casas = zeros(3, 3, 9);
casas(:, :, 1) = [ 1 1 0; 0 0 1; 1 1 0];

casas(:, :, 2) = [ 0 0 0; 0 1 1; 1 1 0];

casas(:, :, 3) = [ 1 1 1; 1 1 0; 1 0 0];

casas(:, :, 4) = [ 1 1 0; 0 0 0; 0 1 1];

casas(:, :, 5) = [ 1 0 1; 0 1 0; 0 0 1];

casas(:, :, 6) = [1 0 1; 1 0 0; 0 0 1];

casas(:, :, 7) = [ 1 0 1; 0 0 0; 1 0 1];

casas(:, :, 8) = [1 0 0; 1 1 0; 0 0 0];

casas(:, :, 9) = [0 0 0; 1 0 0; 0 1 0];
indice_casas = cell(1,length(indice));

for o = 1:length(indice)
  indice_casas {o} = find(casas(:,:,o) ~= 0); % encontra os indices de cada casa infectada dentro da matriz casas de um bairro
endfor
%--------------------------------%


%-----inicializacao aleatoria-----%
function Xijk = inicializarAleatorio(m, J, K, indice)
  Xijk = zeros(m, J, K);
  breeding_grounds = indice(randperm(length(indice))); % gera um vetor de indices aleatorios dos bairros infectados
  k = 1;
  for M = 1:length(breeding_grounds)
    base(M) = randi(m); % selecionar aleatoriamente uma base para cada drone k
  endfor
  for j = 1:length(breeding_grounds)
    Xijk(base(j), breeding_grounds(j), k) = 1;
    k = mod(k, K) + 1; % alterna entre os drones
  endfor
end
%---------------------------------%


%-----funcao objetivo-----%
function [f, f_individual, contador] = calcularObjetivo(Xijk, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas)
  f_individual = zeros(1,K);
  f = 0;
  for k = 1:K
    [base, regiao] = find(Xijk(:,:,k) ~= 0); % encontra a base e a regiao associadas ao drone k
    retorno = false;
    contador(k)=0;
    for j = 1:length(regiao)
      indice_temp = indice_casas{j}; % cria matrizes temporarias com os indices das casas infectadas e a matriz casas (para poder analisar se o drone começa da casa 5 ou nao)
      matriz_temp = casas(:,:,j);
      for o = 1: length(indice_temp)-1
        %fprintf('Drone %d, Regiao %d, Casa %d\n', k, regiao(j), indice_temp(o));
        if j == 1
          f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)); % calculo do tempo gasto para o drone k ir da base i ate o regiao j e despejar o produto na primeira casa
          Dk_bateria(k) = Dk_bateria(k) - ((dij(base(j),regiao(j)) / velocidade(k)) / Dk(k)); % desconta o tempo gasto da bateria restante do drone
          if any(indice_temp == 5) % se a primeira casa for a casa 5 (o ponto de partida do drone eh o centro do bairro j, ou seja, a casa 5)
            f_individual(k) = f_individual(k) + (50 / velocidade(k)) + (ppc / Vk(k)); % calcula o tempo para o drone ir ateh o centro do bairro (casa 5) e dedetiza-lo (50m e a distancia da borda do bairro ateh o centro)
            Dk_bateria(k) = Dk_bateria(k) - (((50 / velocidade(k)) + (ppc / Vk(k))) / Dk(k)); % desconta o tempo gasto da bateria restante do drone
            Tk(k) = Tk(k) - (ppc / Tk(k)); % desconta o produto usado do tanque
            matriz_temp(2,2) = 0; % torna a casa 5 como dedetizada
            while indice_temp(length(indice_temp)) ~= 5 % organiza para que a casa 5 fique no ultimo indice que nao sera acessado pois ja foi acessado
              indice_5 = find(indice_temp == 5); % indica o indice da casa 5
              troca = indice_5 + 1;
              memoria = indice_temp(indice_5);
              indice_temp(indice_5) = indice_temp(troca);
              indice_temp(troca) = memoria;
            endwhile
            indice_temp(length(indice_temp))=0;
          else
            f_individual(k) = f_individual(k) + ((50 + dcc(5,indice_temp(o))) / velocidade(k)) + (ppc / Vk(k)); %calcula o tempo para o drone ir ateh o centro do bairro e em seguida ateh a primeira casa a ser dedetizada e o tempo para dedetiza-la)
            Dk_bateria(k) = Dk_bateria(k) - ((((50 + dcc(5,indice_temp(o)))/ velocidade(k))  + (ppc / Vk(k))) / Dk(k));
            Tk(k) = Tk(k) - (ppc / Tk(k));
          endif
        elseif o > 1 && Tk(k) - (ppc / Tk(k)) >= 0 && Dk_bateria(k) - ((((100 + dcc(indice_temp(o-1), indice_temp(o)) + dij(base(j), regiao(j))) / velocidade(k)) + (ppc / Vk(k))) / Dk(k)) >= 0 % neste caso eh utilizado 100 pois ele precisa retornar ateh o centro para ir ate a borda do bairro novamente e retornar a base
          if retorno == true
            f_individual(k) = f_individual(k) + ((dij(base(j - 1),regiao(j)) + 50 + dcc(5,indice_temp(o))) / velocidade(k)) + (ppc / Vk(k));
            Dk_bateria(k) = Dk_bateria(k) - ((((dij(base(j - 1),regiao(j)) + 50 + dcc(5,indice_temp(o))) / velocidade(k)) + (ppc / Vk(k))) / Dk(k));
            Tk(k) = Tk(k) - (ppc / Tk(k))
            contador(k)=contador(k)+1;
            retorno = false;
          else
            f_individual(k) = f_individual(k) + (dcc(indice_temp(o - 1),indice_temp(o)) / velocidade(k)) + (ppc / Vk(k));
            Dk_bateria(k) = Dk_bateria(k) - (((dcc(indice_temp(o - 1),indice_temp(o)) / velocidade(k)) + (ppc / Vk(k))) / Dk(k));
            Tk(k) = Tk(k) - (ppc / Tk(k));
          endif
        elseif o > 1 && Dk_bateria(k) - ((((100 + dcc(indice_temp(o - 1), indice_temp(o)) + dij(base(j), regiao(j))) / velocidade(k)) + (ppc / Vk(k))) / Dk(k)) <= 0.20 && Tk(k) - (ppc / Tk(k)) <= 0 % se bateria tiver menos que 20% ela eh recarregada (35% da para ir para da base até o local mais longe da cidade, fazer uma casa e voltar)
          f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)) + (tempo_recarga(k) * (1 - Dk_bateria(k))); % calculo do tempo gasto para o drone voltar da regiao atual j ate a base i somado com o tempo de recarga do drone
          Dk_bateria(k) = Dk(k) / Dk(k); % recarrega o drone
          Tk(k) = Tk(k) ./ Tk(k);
          retorno = true;
        elseif o > 1 && Tk(k) - (ppc / Tk(k)) <= 0 && Dk_bateria(k) - ((((100 + dcc(indice_temp(o-1), indice_temp(o)) + dij(base(j), regiao(j))) / velocidade(k)) + (ppc / Vk(k))) / Dk(k)) > 0.20 % se acabou o tanque mas ainda tem mais de 35% da bateria
          f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)) + 180; % se bateria maior que 25% nao recarrega para nao viciar (180 s eh a recarga do tanque)
          Tk(k) = Tk(k) ./ Tk(k);
          retorno = true;
        elseif o > 1 && Dk_bateria(k) - ((((100 + dcc(indice_temp(o - 1), indice_temp(o)) + dij(base(j), regiao(j))) / velocidade(k)) + (ppc / Vk(k))) / Dk(k)) <= 0 % se acabou bateria
          f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)) + (tempo_recarga(k) * (1 - Dk_bateria(k)));
          Dk_bateria(k) = Dk(k) / Dk(k);
          Tk(k) = Tk(k) ./ Tk(k);
          retorno = true;
        endif
      endfor
      if j > 1 && Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dij(base(j),regiao(j)) + 100) / velocidade(k)) / Dk(k)) >= 0 && j < length(regiao)
        f_individual(k) = f_individual(k) + ((djj(regiao(j - 1),regiao(j)) + dj) / velocidade(k));
        Dk_bateria(k) = Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dj) / velocidade(k)) / Dk(k));
      elseif j > 1 && Dk_bateria(k) - (((djj(regiao(j - 1),regiao(j)) + dij(base(j),regiao(j))) / velocidade(k)) / Dk(k)) <=0
        f_individual(k) = f_individual(k) + (dij(base(j),regiao(j)) / velocidade(k)) + (tempo_recarga(k) * (1 - Dk_bateria(k)));
        Dk_bateria(k) = Dk(k) / Dk(k);
        Tk(k) = Tk(k) ./ Tk(k);
        retorno = true;
      elseif j == length(regiao)
        f_individual(k) = f_individual(k) + ((dij(base(j),regiao(j)) + 100) / velocidade(k));
      endif
    endfor
  endfor
  f = sum(f_individual(:));
end
%---------------------------------%


%-----verificar restricoes-----%
function violacao = verificarRestricoes(X_ijk, J, indice)
  violacao = false;
  for j = 1:length(indice)
    total_visitas = sum(sum(X_ijk(:,indice(j),:))); % cada regiao j infectada deve ser visitada por exatamente um drone uma unica vez
    if total_visitas ~= 1
      violacao = true;
      fprintf('Violação na região %d: total de visitas é %d\n', j, total_visitas);
      return;
    endif
  endfor
  if sum(X_ijk(:)) ~= length(indice) % todas as regioes devem ser visitados
    violacao = true;
    return;
  endif
end

% checa se os cada drone eh responsavel por pelo menos 20% da carga de trabalho, ou seja, visita pelo menos 20% dos locais infectados (20% para ter uma margem para o código nao violar as restricoes sempre que tiver um resultado nao perfeito)
function violacao2 = restricao_sobrecarga(X_ijk, J, indice)
  violacao2 = false;
  if sum(sum(X_ijk(:,:,1))) < floor(0.20 * length(indice)) || sum(sum(X_ijk(:,:,2))) < floor(0.20 * length(indice)) || sum(sum(X_ijk(:,:,3))) < floor(0.20 * length(indice)) || sum(sum(X_ijk(:,:,4))) < floor(0.20 * length(indice))
    violacao = true;
    return;
  endif
end

%------------------------------%



%-----crossover (pega parte das regioes do pai 1 com parte das regioes do pai 2)-----%
function filho = crossover(pai1, pai2, m, J, K)
  ponto_corte = randi([1 J-1]); % ponto de corte aleatorio
  filho = zeros(m, J, K); % inicializa o filho
  for k = 1:K
    filho(:,1:ponto_corte,k) = pai1(:,1:ponto_corte,k);
    filho(:,ponto_corte+1:end,k) = pai2(:,ponto_corte+1:end,k);
  endfor
end
%--------------------%


%-----mutacao (troca aleatoria entre 2 bases)-----%
function individuo_mutado = mutacao(individuo, m, J, K, probabilidade_mutacao)
  individuo_mutado = individuo;
  if rand() < probabilidade_mutacao
    % implementacao da mutacao de troca aleatoria de bases
    idx = randi([1 m], 1, 2); % gera dois indices aleatorios para bases
    temp = individuo_mutado(idx(1), :, :); % armazena temporariamente a base escolhida
    individuo_mutado(idx(1), :, :) = individuo_mutado(idx(2), :, :); % troca as bases
    individuo_mutado(idx(2), :, :) = temp; % finaliza a troca
  end
end
%-------------------%


%---selecao cultural---%
function [individuo_aprimorado, f_aprimorado] = buscaLocal(individuo, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas, indice)
  % Copia o indivíduo para aprimorar
  individuo_aprimorado = individuo;

  for k = 1:K
    % Encontra as regiões visitadas pelo drone k
    [base, regiao] = find(individuo(:,:,k) ~= 0);

    if length(regiao) > 1
      % Seleciona aleatoriamente duas regiões para trocar
      i = randi(length(regiao));
      j = randi(length(regiao));
      while j == i
        j = randi(length(regiao));
      endwhile

      % Troca as regiões
      temp = regiao(i);
      regiao(i) = regiao(j);
      regiao(j) = temp;

      % Atualiza o indivíduo com a nova ordem
      individuo_aprimorado(:,:,k) = 0;
      for idx = 1:length(regiao)
        individuo_aprimorado(base(idx), regiao(idx), k) = 1;
      endfor
    endif
  endfor

  % Recalcula a função objetivo para garantir que a nova solução é válida
  if verificarRestricoes(individuo_aprimorado, J, indice)
    [f_aprimorado, ~] = Inf;
    [f_original, ~] = calcularObjetivo(individuo, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
  elseif restricao_sobrecarga(individuo_aprimorado, J, indice)
    [f_aprimorado, ~] = calcularObjetivo(individuo_aprimorado, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    [f_original, ~] = calcularObjetivo(individuo, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    f_aprimorado += 600;
  else
    [f_aprimorado, ~] = calcularObjetivo(individuo_aprimorado, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    [f_original, ~] = calcularObjetivo(individuo, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    % Substitui o indivíduo original apenas se a nova solução for melhor
    if f_aprimorado < f_original
      individuo_aprimorado = individuo_aprimorado;
    else
      individuo_aprimorado = individuo;
    endif
  endif
end
%----------------------%


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
tic;
probabilidade_mutacao = 0.05;

% Inicializa a população
populacao = zeros(m, J, K, pop_size);
f = zeros(pop_size, 1);

for i = 1:pop_size
  populacao(:,:,:,i) = inicializarAleatorio(m, J, K, indice);
endfor

% Loop de gerações
for iter = 1:geracoes
  Dk_bateria = Dk ./ Dk;
  Tk_tanque = Tk ./ Tk;
  % Avalia a população
  for i = 1:pop_size
    % Verifica restrições
    if verificarRestricoes(populacao(:,:,:,i), J, indice)
      f(i) = Inf; % penaliza solucoes que violam as restricoes
    elseif restricao_sobrecarga(populacao(:,:,:,i), J, indice)
      [f(i), ~] = calcularObjetivo(populacao(:,:,:,i), dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
      f(i) += 600;
    else
      [f(i), ~] = calcularObjetivo(populacao(:,:,:,i), dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    endif
  endfor
  if iter == 1
    inicial = 1
    min(f)
    best_fitness(iter) = min(f);
    atual = min(f);
    Y(iter) = best_fitness(iter);
  elseif min(f) < atual
    iter
    printf("%.4f\n", min(f))
    best_fitness(iter + 1) = min(f);
    atual = minf(f);
  else
    best_fitness(iter + 1) = best_fitness(iter);
  endif

  % atualizacao da chamada da funcao de selecao
  num_pares = pop_size / 2; % Define o número de pares de pais (5 neste caso)
  pais = selecao(populacao, f, num_pares);

  % crossover e mutaçcao
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

    % aplicar mutacao se chance for cumprida
    if probabilidade_mutacao > rand
      filho1 = mutacao(filho1, m, J, K, probabilidade_mutacao);
      filho2 = mutacao(filho2, m, J, K, probabilidade_mutacao);
    endif

      % avalia novos individuos
    if verificarRestricoes(filho1, J, indice)
      f1 = Inf; % penaliza solucoes que violam as restricoes
    elseif restricao_sobrecarga(filho1, J, indice)
      [f1, ~] = calcularObjetivo(filho1, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
      f1 += 600;
    else
      [f1, ~] = calcularObjetivo(filho1, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    endif

    if verificarRestricoes(filho2, J, indice)
      f2 = Inf; % penaliza solucoes que violam as restricoes
    elseif restricao_sobrecarga(filho2, J, indice)
      [f2, ~] = calcularObjetivo(filho2, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
      f2 += 600;
    else
      [f2, ~] = calcularObjetivo(filho2, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    endif
    if f1 < atual && f1 < f2
      iter
      printf("%.4f\n", f1)
      best_fitness(iter + 1) = f1;
      atual = f1;
    elseif f2 < atual && f2 < f1
      iter
      printf("%.4f\n", f2)
      best_fitness(iter + 1) = f2;
      atual = f2;
    endif

    % aplicar busca local aos filhos
    sobrecarga = true;
    while sobrecarga == true
    sobrecarga = false;
    count = count + 1;
    [filho_busca1, f_busca1] = buscaLocal(filho1, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas, indice);
    [filho_busca2, f_busca2] = buscaLocal(filho2, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas, indice);
    for r = 1:K
      if (sum(sum(filho1(:,:,r))) < floor(0.20 * J) && count < 10) || (sum(sum(filho2(:,:,r))) < floor(0.20 * J) && count < 10)
        sobrecarga = true;
      elseif count == 10
        sobrecarga = false;
      endif
    endfor
  endwhile

  % avalia filhos apos selecao cultural
  if f_busca1 < atual && f_busca1 < f_busca2
    iter
    printf("%.4f\n", f_busca1)
    best_fitness(iter + 1) = f_busca1;
    atual = f_busca1;
  elseif f_busca2 < atual && f_busca2 < f_busca1
    iter
    printf("%.4f\n", f_busca2)
    best_fitness(iter + 1) = f_busca2;
    atual = f_busca2;
  endif

  % seleciona os melhores filhos para a proxima populacao (entre individuo cruzado e individuo sob selecao cultural)
  if f_busca1 < f1
    filho1 = filho_busca1;
  endif
    if f_busca2 < f2
    filho2 = filho_busca2;
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
title('Curva de convergência AM');
xlabel('Geração');
ylabel('Valor obtido');
