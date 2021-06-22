% Interpoliavimas_Ermito_splainais
% kiekvienas intervalas tarp tasku interpoliuojamas
% 2 eiles Ermito daugianariais, nustatant isvestiniu reiksmes pagal 
% skaitinio diferencijavimo formules

function Ermito_splainai_pagal_nutylejima
clc,close all

% X=[0  0.2  3  4  9  10.5  14 14.5]
% Y=[0.2 1 1.2 3  3  1.2  1  0.2]
X=[0 0.5 2 3 4 5.5 6 7 8]
% Y=[0 1 4 16 25 16 4 1 0]
Y=5./(1+4*(X-4).^2)
% DY=[0  0  0  0  0  0  0  0]
DY=Akima(X,Y)  % Isvestiniu reiksmiu interpoliavimo taskuose nustatymas pagal Akima formules
nP=length(X) % interpoliavimo tasku skaicius

for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
    nnn=100;
    xxx=[X(iii):(X(iii+1)-X(iii))/nnn:X(iii+1)];
    fff=0;
    for j=1:2
        [U,V]=Hermite(X(iii:iii+1),j,xxx);
        fff=fff+U*Y(iii+j-1)+V*DY(iii+j-1);
    end
    figure(1), hold on, grid on, axis equal
    plot(xxx,fff,'r-','LineWidth',2.5);
    plot(X(iii:iii+1),Y(iii:iii+1),'ko','LineWidth',2,'MarkerSize',8)
end %-----------------ciklo per intervalus pabaiga
legend({sprintf('Ermito splainai %d intervaluose',nP-1)});

return
end

function [U,V]=Hermite(X,j,x)   % Ermito daugianariai
    L=Lagrange(X,j,x); DL=D_Lagrange(X,j,X(j));
    U=(1-2*DL.*(x-X(j))).*L.^2;
    V=(x-X(j)).*L.^2;
return
end

function L=Lagrange(X,j,x) % Lagranzo daugianaris
    n=length(X);
    L=1;
    for k=1:n, if k ~= j, L=L.*(x-X(k))/(X(j)-X(k)); end,end
return
end

function DL=D_Lagrange(X,j,x) % Lagranzo daugianario isvestine pagal x
    n=length(X);
    DL=0; %DL israskos skaitiklis
    for i=1:n % ciklas per atmetamus narius
        if i==j, continue, end 
        Lds=1;
        for k=1:n, if k ~= j && k ~= i , Lds=Lds.*(x-X(k)); end,  end
        DL=DL+Lds;
    end
    Ldv=1;   %DL israskos vardiklis 
    for k=1:n, if k ~= j, Ldv=Ldv.*(X(j)-X(k)); end, end
    DL=DL/Ldv;
return
end

function DY=Akima(X,Y) % Isvestiniu reiksmiu interpoliavimo taskuose nustatymas pagal skaitinio integravimo formules
    n=length(X);
    f=inline('(2*x-xi-xip1)/((xim1-xi)*(xim1-xip1))*yim1+(2*x-xim1-xip1)/((xi-xim1)*(xi-xip1))*yi+(2*x-xim1-xi)/((xip1-xim1)*(xip1-xi))*yip1')
    for i=1:n
        if i == 1,xim1=X(1);xi=X(2);xip1=X(3); yim1=Y(1);yi=Y(2);yip1=Y(3);DY(i)=f(xim1,xi,xim1,xip1,yi,yim1,yip1);
        elseif i == n, xim1=X(n-2);xi=X(n-1);xip1=X(n); yim1=Y(n-2);yi=Y(n-1);yip1=Y(n); DY(n)=f(xip1,xi,xim1,xip1,yi,yim1,yip1);
        else, xim1=X(i-1);xi=X(i);xip1=X(i+1); yim1=Y(i-1);yi=Y(i);yip1=Y(i+1); DY(i)=f(xi,xi,xim1,xip1,yi,yim1,yip1);
        end
    end
return
end
