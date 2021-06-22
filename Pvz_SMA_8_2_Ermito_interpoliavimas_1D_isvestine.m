% Ermito_interpoliavimas
% Interpoliuojancio daugianario isvestine

function Ermito_interpoliavimas_isvestine
clc,close all
syms  f x 

f=sin(2*x), df=diff(f)

nP=5
xrange=[-pi,pi]
nnn=3000
X=[xrange(1):(xrange(2)-xrange(1))/(nP-1):xrange(2)] 
xxx=[xrange(1):(xrange(2)-xrange(1))/nnn:xrange(2)];

x=X;
F=eval(f)
DF=eval(df)

n=length(X)

% return
spalvos=['b','r','g','c','m','k','b','r','g','c','m','k'];
figure(1); hold on; grid on; title('Ermito daugianariai U');axis equal
figure(2); hold on; grid on; title('Ermito daugianariai V');axis equal
figure(3); hold on; grid on; title('Ermito daugianariu isvestines U');axis equal
figure(4); hold on; grid on; title('Ermito daugianariu isvestines V');axis equal

fff=0;dfff=0;
legU={};legV={};
for j=1:n
    [U,V]=Hermite(X,j,xxx);
    [DU,DV]=D_Hermite(X,j,xxx);
    figure(1);plot(xxx,U,[spalvos(j),'-']);
    figure(2);plot(xxx,V,[spalvos(j),'-']);
    figure(3);plot(xxx,DU,[spalvos(j),'-']);
    figure(4);plot(xxx,DV,[spalvos(j),'-']);
    legU{j}=sprintf('U%2d',j);
    legV{j}=sprintf('V%2d',j);
    legDU{j}=sprintf('DU%2d',j);
    legDV{j}=sprintf('DV%2d',j);
    fff=fff+U*F(j)+V*DF(j);
    dfff=dfff+DU*F(j)+DV*DF(j);
end

figure(1);legend(legU);plot(X,zeros(1,n),'k*')
figure(2);legend(legV);plot(X,zeros(1,n),'k*')
figure(3);legend(legDU);plot(X,zeros(1,n),'k*')
figure(4);legend(legDV);plot(X,zeros(1,n),'k*')


figure(5), hold    on, grid on, axis equal
x=xxx;plot(xxx,eval(f),'b-');plot(xxx,fff,'c-','LineWidth',2.5);
plot(xxx,eval(df),'r-');plot(xxx,dfff,'m-','LineWidth',2.5);

legend({'duotoji funkcija sin(2x)',[' Ermito interpoliavimas per ',sprintf('%d',n),' taskus f'],...
    'duotosios funkcijos isvestine',[' Ermito interpoliavimas per ',sprintf('%d',n),' taskus df/dx']})
plot(X,F,'ko','LineWidth',2,'MarkerSize',8)

return
end

function [U,V]=Hermite(X,j,x)
% Ermito daugianariai
n=length(X);
L=Lagrange(X,j,x); DL=D_Lagrange(X,j,X(j));
U=(1-2*DL.*(x-X(j))).*L.^2;
V=(x-X(j)).*L.^2;

return
end

function [DU,DV]=D_Hermite(X,j,x)
% Ermito daugianariu isvestines
n=length(X);
L=Lagrange(X,j,x); 
DLj=D_Lagrange(X,j,X(j));
DL=D_Lagrange(X,j,x);

DU=-2*DLj*L.^2+2*(1-2*DLj*(x-X(j))).*L.*DL;
DV=L.^2+2*(x-X(j)).*L.*DL;

return
end


function L=Lagrange(X,j,x)
n=length(X);
% Lagranzo daugianaris
L=1;
for k=1:n
    if k ~= j, L=L.*(x-X(k))/(X(j)-X(k)); end
end

return
end


function DL=D_Lagrange(X,j,x)
% Lagranzo daugianario isvestine pagal x
n=length(X);
DL=0; %DL israskos skaitiklis
for i=1:n % ciklas per atmetamus narius
    if i==j, continue, end 
    Lds=1;
    for k=1:n  
        if k ~= j && k ~= i , Lds=Lds.*(x-X(k)); end
    end
    DL=DL+Lds;
end
    Ldv=1;   %DL israskos vardiklis 
    for k=1:n 
        if k ~= j, Ldv=Ldv.*(X(j)-X(k)); end 
    end
DL=DL/Ldv;

return
end
