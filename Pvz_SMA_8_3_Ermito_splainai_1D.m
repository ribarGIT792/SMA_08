% Interpoliavimas_Ermito_splainais
% kiekvienas intervalas tarp tasku interpoliuojamas
% 2 eiles Ermito daugianariais

function Ermito_splainai
clc,close all
syms  f x 

% f=sin(x)  % duotoji funkcija
f=1./(1+5*x.^2)
df=diff(f)
nP=5 % interpoliavimo tasku skaicius
xrange=[-pi,pi]
X=[xrange(1):(xrange(2)-xrange(1))/(nP-1):xrange(2)] 
F=eval(subs(f,sym(x),sym(X)))
DF=eval(subs(df,sym(x),sym(X)))

for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
    nnn=100;
    xxx=[X(iii):(X(iii+1)-X(iii))/nnn:X(iii+1)];
    fff=0;
    for j=1:2
        [U,V]=Hermite(X(iii:iii+1),j,xxx);
        fff=fff+U*F(iii+j-1)+V*DF(iii+j-1);
    end
    figure(1), hold on, grid on, axis equal
    plot(xxx,eval(subs(f,sym(x),sym(xxx))),'b-');
    plot(xxx,fff,'r-','LineWidth',2.5);
    plot(X(iii:iii+1),F(iii:iii+1),'ko','LineWidth',2,'MarkerSize',8)
end %-----------------ciklas per intervalus pabaiga
legend({['duotoji funkcija',char(f)],sprintf('Ermito splainai %d intervaluose',nP-1)});
return
end

function [U,V]=Hermite(X,j,x)  % Ermito daugianariai
    L=Lagrange(X,j,x); DL=D_Lagrange(X,j,X(j));
    U=(1-2*DL.*(x-X(j))).*L.^2;
    V=(x-X(j)).*L.^2;
return
end

function L=Lagrange(X,j,x)   % Lagranzo daugianaris
    n=length(X);
    L=1;
    for k=1:n, if k ~= j, L=L.*(x-X(k))/(X(j)-X(k)); end, end
return
end

function DL=D_Lagrange(X,j,x)  % Lagranzo daugianario isvestine pagal x
    n=length(X);
    DL=0; %DL israskos skaitiklis
    for i=1:n % ciklas per atmetamus narius
        if i==j, continue, end 
        Lds=1;
        for k=1:n, if k ~= j && k ~= i , Lds=Lds.*(x-X(k)); end, end
        DL=DL+Lds;
    end
    Ldv=1;   %DL israskos vardiklis 
    for k=1:n, if k ~= j, Ldv=Ldv.*(X(j)-X(k)); end,  end
    DL=DL/Ldv;
return
end
