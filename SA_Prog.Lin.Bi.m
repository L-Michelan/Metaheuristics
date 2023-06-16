clear all
close all
clc
tic
n=4; %tamanho do vetor solucao%
l=101;
L=101 %criterio de parada kmax%
k=1
v=rand;
%-----Geração do vetor soluçao inicial-----%
for i=1:n
  Q(i)=i;
  r=rand;
  if r>0.5
    x(i)=1;
  else
    x(i)=0;
  endif
  E(i)=x(i);%vetor elite%
endfor
%--------------------Fim-------------------%

%---avaliacao da solucaoao inicial, maximizar f---%
f=8*x(1)+5*x(2)+6*x(3)+4*x(4);
if x(3)+x(4)>1
  F=0;
elseif x(4)>x(2)
  F=0;
  elseif x(3)>x(4)
    F=0;
    elseif 6*x(1)+3*x(2)+5*x(3)+2*x(4)>10
      F=0;
      else
      F=1;
endif
c=f+50*F;
%----------------------Fim------------------------%
H=0;
%---contador(explora toda vizinhanca da solucao aleatoriamente até encontrar uma solucao melhor)---%
while k<=l
  Q(H+1)=H+1;
  H=H+1;
%--------Manipulador de vizinhanca--------%
p=randi(4); %operador manipulador da vizinhanca%
P=randi(5);
h=P-p;
if h<0
  h=p-P;
elseif h==0
  h=1;
  endif
A=x(p);
x(p)=x(h);
x(h)=A;
%-------------------Fim-------------------%


%---avaliacao da nova solucaoao---%
j=8*x(1)+5*x(2)+6*x(3)+4*x(4);
if x(3)+x(4)>1
  F=0;
elseif x(4)>x(2)
  F=0;
  elseif x(3)>x(4)
    F=0;
    elseif 6*x(1)+3*x(2)+5*x(3)+2*x(4)>10
      F=0;
      else
      F=1;
endif
d=j+50*F
if (d==63)
  k=l;
  V(H)=d;
  endif
if d>c
  for i=1:n
  E(i)=x(i);
  V(H)=d;
endfor
elseif v=L/400
for i=1:n
E(i)=x(i);
V(H)=d;
endfor
endif
%--------------Fim----------------%
k=k+1;
L=L-1;
endwhile
%----------------------------------------------Fim--------------------------------------------------%
Time=toc;
plot(Q,V)
xlabel('GERAÇOES')
ylabel('VALOR OBTIDO')
title('gráfico do numero de geraçoes pelo valor obtido')
xlswrite('V1.xlsx',V)
xlswrite('e1.xlsx',E)
xlswrite('Time.xlsx',Time)
