clear all
clc

%definir uma matriz fixa para bairros e casas

%-----Parametros e variaveis-----%
max_iter = 100; % número máximo de iterações
Y = zeros(101);% variavel para o grafico e resultado
X = 1:1:max_iter+1;

J = 16; % total de locais da cidade a serem visitados
K = 4; % numero de drones
m = 8; % numero de locais potenciais para a instalacao das bases

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

dij = dij'; %toda matriz que estao transposta serve apenas para que os indices da matriz referencie exatamente aos indices da imagem da cidade hipotetica
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
%-------------------------%


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
%------------------------------%


%-----busca local (troca de bases entre dois drones)-----%
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


%-----mudanca de vizinhanca (troca as matrizes entre os drones)-----%
function Xijk = mudarVizinhanca(Xijk, K)
  % seleciona dois drones aleatorios
  drone1 = randi(K);
  drone2 = randi(K);
  while drone2 == drone1
    drone2 = randi(K); % garante que os drones selecionados sejam diferentes
  endwhile

  % troca todas as bases entre os drones mantendo as regioes inalteradas
  for i = 1:size(Xijk, 1)
    % salva as visitas das bases atuais
    visitas_base_d1 = Xijk(i, :, drone1);
    visitas_base_d2 = Xijk(i, :, drone2);

    % troca as visitas das bases entre os drones
    Xijk(i, :, drone1) = visitas_base_d2;
    Xijk(i, :, drone2) = visitas_base_d1;
  endfor
end
%---------------------------------------------------------------%


%-----Algoritmo VNS-----%
function [X_best, f_best, f_individual_best, contador, Y] = VNS(dj, Dk, m, J, K, max_iter, dij, djj, tempo_recarga, velocidade, Tk, Vk, Vp, dcc, ppc, casas, indice_casas, indice)

  % inicializacaoo aleatoria
  Dk_bateria = Dk ./ Dk; % inicia a bateria (iniciada como 1 para representar 100%)
  Tk_tanque = Tk ./ Tk; % inicia o tanque
  iniciar=0;
  Xijk = inicializarAleatorio(m, J, K, indice);
  while verificarRestricoes(Xijk, J, indice)
    Xijk = inicializarAleatorio(m, J, K, indice);
  endwhile
  inicial = 1
  [f_best, f_individual_best, contador] = calcularObjetivo(Xijk, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas) % calcula o valor de tempo da primeira solucao gerada para todos os drones e cada drone indiviudual
  X_best = Xijk; % atribui um valor inicial a X_best
  f_curr = Inf;
  count=contador;
  Y(1) = f_best;

  % VNS loop
  for iter = 2:max_iter + 1
    Dk_bateria = Dk ./ Dk; % reseta a bateria dos drones
    Tk_tanque = Tk ./ Tk; % reseta o tanque

    % inicializacao da solucao
    X_ijk = X_best;

    % inicializa uma celula para armazenar as solucoes da busca local
    X_ijk_temp = cell(1, 3);
    % executa a busca local 3 vezes e seleciona a melhor solucao perturbada
    equal = 0;
    for o = 1:3
      X_ijk_temp{o} = buscaLocal(X_ijk, K);

      % verifica se a nova solucao e igual a alguma das solucoes anteriores
      if o ~= 1
        for i = 1:o-1
          while isequal(X_ijk_temp{o}, X_ijk_temp{i}) && equal <=20
            equal = equal + 1;
            X_ijk_temp{o} = buscaLocal(X_ijk, K); % Gera nova solução se for igual a alguma anterior
          endwhile
        endfor
      endif

      % calcula o objetivo para a nova solucao
      [f_temp(o), f_individual_temp(o,:), contador] = calcularObjetivo(X_ijk_temp{o}, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);

      % atualiza a melhor solucao se a nova for melhor
      if o == 1 || f_temp(o) < f_curr
        f_curr = f_temp(o);
        f_individual_curr = f_individual_temp(o,:);
        X_ijk = X_ijk_temp{o};
        count= contador;
      endif
    endfor

    % compara solucao perturbada com a melhor solucao
    if f_curr < f_best && ~verificarRestricoes(X_ijk, J, indice)
      iter
      f_curr
      X_best = X_ijk;
      f_best = f_curr;
      f_individual_best = f_individual_curr;
      count=contador;
    endif

    % mudança de vizinhancaa
    X_ijk = mudarVizinhanca(X_ijk, K);

    % avaliacao da solucao apos mudanca de vizinhanca
    [f_curr, f_individual_curr, contador] = calcularObjetivo(X_ijk, dj, K, dij, djj, Dk_bateria, Dk, tempo_recarga, velocidade, Tk, Tk_tanque, Vk, Vp, dcc, ppc, casas, indice_casas);
    if f_curr < f_best && ~verificarRestricoes(X_ijk, J, indice)
      iter
      f_curr
      X_best = X_ijk;
      f_best = f_curr;
      f_individual_best = f_individual_curr;
      count=contador;
    endif
    Y(iter) = f_best;
  endfor
end
%-----------------------%


% chamar o algoritmo VNS
tic;
[X_best, f_best, f_individual_best, count, Y] = VNS(dj, Dk, m, J, K, max_iter, dij, djj, tempo_recarga, velocidade, Tk, Vk, Vp, dcc, ppc, casas, indice_casas, indice);
tempo_CPU = toc

% exibir resultado
disp('Melhor solução encontrada:');
disp(X_best);
disp(['Tempo total: ', num2str(f_best)]);
disp(['Tempo de cada drone: ', num2str(f_individual_best)]);

plot(X,Y);
title('Curva de convergência VNS');
xlabel('Iteração');
ylabel('Valor obtido');
