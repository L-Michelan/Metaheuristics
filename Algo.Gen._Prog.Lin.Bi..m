clear all
close all
clc
tic
f=0; %vetor fitness%
n=4 %input ('Tamanho do vetor ')%
p=6 %input ('Tamanho da popula��o ')%
a=0; %matriz popula�ao%
c=0;
R=0; %vetor restri�ao 0=nao 1=sim
s=1+round(rand*p) %tamanho popula�ao para cruzamento%
b=0; %matriz popula�ao cruzamento%
P=0; %popula��o aumentada
H=[0,0,0,0;0,0,0,0];
E=0;
r=rand;
F=p+4;
V=p+2;
X=100 %input (' N�mero de gera��es ')
r=rand;
U=0
cont=1;


%---transforma variaveis em matrizes---%
for i=1:p
  for j=1:n
    a(i,j)=0;
    B(j)=0;
    E(i,j)=0;
    e(i)=0;
  j=j+1;
endfor
c(i)=0;
R(i)=0;
f(i)=0;
i=i+1;
endfor
%-----------------Fim------------------%

%---gera matriz pop inicial---%
for i=1:p
  for j=1:n
    h=rand
  if h>0.5
    a(i,j)=1
  else
    a(i,j)=0
  endif
  j=j+1
endfor
i=i+1
endfor
%------------Fim-------------%
A=1;
%---criterio de parada---%
while A<X %crit�rio de parada
r=rand;
  %---gera vetor pop aumentada---%
if r<0.05
for i=1:F
  for j=1:n
    P(i,j)=0;
  endfor
endfor
else
for i=1:p+2
  for j=1:n
    P(i,j)=0;
  endfor
endfor
endif
%--------------Fim--------------%

  U

%---avaliacao pop inicial---%
for i=1:p
c(i)=8*a(i,1)+5*a(i,2)+6*a(i,3)+4*a(i,4);#func obj
if a(i,3)+a(i,4)>1
  R(i)=0;
elseif a(i,3)>a(i,4)
  R(i)=0;
  elseif a(i,4)>a(i,2)
    R(i)=0;
    elseif 6*a(i,1)+3*a(i,2)+5*a(i,3)+2*a(i,4)> 10
      R(i)=0;
     else
       R(i)=1;
       endif
  f(i)=c(i)+50*R(i);
  if f(i)==63
  A=X+1;
  endif
if f(i)>e(i) #avalia elite
  for j=1:n
  e(i)=f(i);
  E(i,j)=a(i,j);
  endfor
  endif
endfor
%------------Fim------------%

T=0;

%---cria matriz com novo tamanho de populacao---%
if (s==0)
  s=1
  endif
for l=1:s
  t=round(rand*p);
  if (t>=p)
    t=t-2;
  elseif (t==0)
    t=t+1;
  endif
  for o=1:n
    b(l,o)=a(t,o);
endfor
endfor
%---------------------Fim-----------------------%

%---cruzamento---%
for o=1:2 #seleciona  ind da nova pop p/ cruzar
  T=round(rand(1)*s);
  if T==0
    T=T+1;
  endif
  for Z=1:n
    B(o,Z)=b(T,Z);
  endfor
endfor
t=4;
t=round(rand*n);
T=1;
if (t>=n) %cruzamento%
    t=t-1;
  elseif (t<=1)
    t=2;
  endif
  if (s==1)
    I=B(1,t+1);
    B(1,t+1)=B(1,t-1);
    B(1,t-1)=I;
else
I=B(1,t+1);
B(1,t+1)=B(2,t+1);
B(2,t+1)=I;
T=t
endif
if t>1
I=B(1,t-1);
B(1,t-1)=B(2,t-1);
B(2,t-1)=I;
elseif t<=1
I=B(1,t+2);
B(1,t+2)=B(2,t+2);
B(2,t+2)=I;
endif
%------Fim-------%

%---avaliacao prole---%
for i=1:2
c(i)=8*B(i,1)+5*B(i,2)+6*B(i,3)+4*B(i,4);
if B(i,3)+B(i,4)>1
  R(i)=0;
elseif B(i,3)>B(i,4)
  R(i)=0;
  elseif B(i,4)>B(i,2)
    R(i)=0;
    elseif 6*B(i,1)+3*B(i,2)+5*B(i,3)+2*B(i,4)> 10
      R(i)=0;
     else
       R(i)=1;
      endif
      for o=1:p;
f(i)=c(i)+50*R(i);
if f(i)==63
  A=X+1;
  endif
if f(i)>e(o)
    for j=1:n
  e(o)=f(i);
  E(o,j)=B(i,j);
endfor
endif
endfor
endfor
%-------------Fim-------------%

r

%---mutacao---%
if r<0.05%muta��o%
  for L=1:2
    v=rand;
    for j=1:n
  if (v>0.5)
    H(L,j)=B(L,j);
    B(L,j)=1;
else
  H(L,j)=B(L,j);
  B(L,j)=0;
  endif
    endfor
  endfor
%-----Fim-----%

%---avaliacao idiv mutados---%
for i=1:2
c(i)=8*B(i,1)+5*B(i,2)+6*B(i,3)+4*B(i,4);%vetor fun�ao obj%
if B(i,3)+B(i,4)>1
  R(i)=0;
elseif B(i,3)>B(i,4)
  R(i)=0;
  elseif B(i,4)>B(i,2)
    R(i)=0;
    elseif 6*B(i,1)+3*B(i,2)+5*B(i,3)+2*B(i,4)> 10
      R(i)=0;
     else
       R(i)=1;
      endif
      for o=1:p
f(i)=c(i)+50*R(i);
if f(i)==63
  A=X+1;
  endif
if f(i)>e(o)
    for j=1:n
  e(o)=f(i);
  E(o,j)=B(i,j);
    endfor
endif
endfor
endfor
endif
%------------Fim-------------%

%---gera pop gerada U pop indiv U pop mutacao, se mutou---%
for i=1:p #pop aumentada
  for j=1:n
    P(i,j)=a(i,j);
  endfor
endfor
for i=p:p+2 #pop aumentada 2
  for l=1:2
    for j=1:n
    P(i,j)=B(l,j);
  endfor
endfor
endfor
if (r<0.05) #pop aumentada se mutou
  for i=p+2:F
    for l=1:2
      for j=1:n
      P(i,j)=H(l,j);
    endfor
  endfor
endfor
%---------------------------Fim----------------------------%

%---avaliacao pop aumentada mutando---%
for i=1:F
c(i)=8*P(i,1)+5*P(i,2)+6*P(i,3)+4*P(i,4);
if P(i,3)+P(i,4)>1
  R(i)=0;
elseif P(i,3)>P(i,4)
  R(i)=0;
  elseif P(i,4)>P(i,2)
    R(i)=0;
    elseif 6*P(i,1)+3*P(i,2)+5*P(i,3)+2*P(i,4)> 10
      R(i)=0;
     else
       R(i)=1;
      endif
      for L=1:p
f(i)=c(i)+50*R(i);
if f(i)==63
  A=X+1;
endif
if f(i)>e(L)
    for j=1:n
  e(L)=f(i);
  E(L,j)=P(i,j);
    endfor
endif
endfor
endfor
%-----------------Fim------------------%

%---avaliacao pop aumentada nao mutando---%
elseif r>=0.05
X
for i=1:V
  c(i)=8*P(i,1)+5*P(i,2)+6*P(i,3)+4*P(i,4);
if P(i,3)+P(i,4)>1
  R(i)=0;
elseif P(i,3)>P(i,4)
  R(i)=0;
  elseif P(i,4)>P(i,2)
    R(i)=0;
    elseif 6*P(i,1)+3*P(i,2)+5*P(i,3)+2*P(i,4)>10
      R(i)=0;
     else
       R(i)=1;
      endif
      for o=1:p
f(i)=c(i)+50*R(i);
if f(i)==63
  A=X+1;
  endif
if f(i)>e(o)
  for j=1:n
  e(o)=f(i);
  E(o,j)=P(i,j);
  endfor
endif
endfor
endfor
endif
%-----------------Fim---------------------%

%---nova geracao de pop---%
for i=1:p
  t=1+round(rand*F-3);
  if t>=F
    t=F-2;
  elseif t<=0
    t=1;
    endif
  for j=1:n
    a(i,j)=P(t,j);
  endfor
endfor
%-----------Fim-----------%

U=U+1;
a;
B;
E;
s=1+round(rand*p);
F=p+4;
V=p+2;
r=rand;

%---reseta pop aumentada se mutou ou n---%
 if r<0.05
for I=1:F
  for J=1:n
    P(I,J)=0;
  endfor
endfor
else
for I=1:V
  for J=1:n
    P(I,J)=0;
  endfor
endfor
endif
%------------------Fim-------------------%
Q(cont)=cont;
e
elite(cont)=max(e);
cont=cont+1;
A=A+1;
endwhile
Time=toc;
%----------Fim-----------%
plot(Q,elite)
xlabel('GERAÇOES')
ylabel('ELITE DA GERAÇAO')
title('gráfico do numero de geraçoes pelo valor obtido')
xlswrite('avaliado.xlsx',elite)
xlswrite('Pop.xlsx',E)
xlswrite('time.xlsx',Time)
