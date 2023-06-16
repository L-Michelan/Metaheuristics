#f=2*exp(x)-x*sin(x+3)-3// função a ser aplicada a busca de zeros para todo -120<x<0
%variaveis ja usadas: a_A_b_B_c_C_d_D_e_E_f_F_G_t_T_i_j_J_k_K_p_n_R
#plotar grafico de elite obtida por iteracao
clear all
close all
clc
a=0
b=0;
c=0;
k=100;
d=0;
i=0;
j=0;
E=0;
K=0;
d=1;
C=0#contador de zeros
Q=[1:1:k+1];
for i=1:37
E(i)=1; #vetor específico para os zeros da funcao
endfor
for i=1:k+1
e(i)=-1; #vetor elite da prole
endfor
%---gerar vetores solução, onde cada a(i) representa um individuo---%
for i=1:10
  a(i)=-rand*120
endfor
for i=1:10
  f=2*exp(a(i))-a(i)*sin(a(i)+3)-3
  if (f>0)
  A(i)=f; %vetor resultado
  elseif (f<0)
  A(i)=-f;
  else
    E(1)=a(i);
    A(i)=f;
    C=C+1;
    d=d+1;
    endif
endfor
e(1)=-min(A);
%-------------------------------------------------------------------%

%---organiza indiv em ordem crescente para separar melhores depois---%
n=length(a);
m=true;
while m==true
        m=false
        for i=1:(n-1)
            if A(i) > A(i+1)
                aux=A(i);
                aux2=a(i);
                A(i)=A(i+1);
                a(i)=a(i+1);
                A(i+1)=aux;
                a(i+1)=aux2;
                m=true;
            endif
        endfor
    endwhile
    R(1)=a(1); #R é o vetor q indica o individuo que possuiu o melhor fitness da geracao
%-------------------------------------------------------------%

%---Critério de parada---%
for K=1:k
%---cruzamento---%

  for i=1:10
    dominio=true;
    while dominio==true
      dominio=false;
r=rand;
h=randi(10);
g=1+rand*2
    if (r>0.5)
    c(i)=a(h)*g; #trabalhando com multiplicacao e divisao para evitar sair do domínio [0, -120]
  else
    c(i)=a(h)/g;
  endif
  if (c(i)<-120)
    dominio=true
  endif
  endwhile
%----------------%

%---avaliação da prole---%
    F=2*exp(c(i))-c(i)*sin(c(i)+3)-3;
    if (F>0)
      p(i)=F;
    elseif (F<0)
      p(i)=-F;
    else
      p(i)=0;
      E(K+1)=c(i);
      C=C+1;
    endif
    e(K+1)=-min(p);
  endfor
%------------------------%

%---organiza indiv em ordem crescente para separar melhores depoios---%
n=length(c);
m=true;
while m==true
        m=false
        for i=1:(n-1)
            if p(i)>p(i+1)
                aux3=p(i);
                aux4=c(i);
                p(i)=p(i+1);
                c(i)=c(i+1);
                p(i+1)=aux3;
                c(i+1)=aux4
                m=true
            endif
        endfor
    endwhile
    R(K+1)=c(1)
%----------------------------------------------%

%---cria populacao com os 5 melhores das pop 1 e 2---%
for i=1:10
  if (i<=5)
  D(i)=c(i);
elseif (i>5)
  for j=1:5
  D(i)=a(j);
endfor
endif
  endfor
%-------------------------------------------------------%

%---analisa para ver se ficou algum outro melhor para tras---%
for i=1:10
for j=5:10
if (i<=5)
  if (D(i)<=A(j))
    D(i)=a(j);
  endif
else
  if (D(i)<=p(j))
    D(i)=c(j);
  endif
endif
endfor
endfor
%------------------------------------------------------------%

%---substitui a pop inicial pela pop dos mais aptos e reseta pop mais aptos---%
for i=1:10
  a(i)=D(i)
  D(i)=0;
  endfor
%-----------------------------------------------------------------------------%

endfor
%------------------------%
xlswrite('R1.xlsx',R)
xlswrite('e1.xlsx',e)
xlswrite('Q1.xlsx',Q)
plot3(Q,e,R)
xlabel('GERAÇOES')
ylabel('VALOR OBTIDO')
title('gráfico do numero de geraçoes pelo valor obtido')
