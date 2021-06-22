 Ermito_interpoliavimas

function Ermito_interpoliavimas
clc,close all
syms  f x 

f=sin(x)
f =1./(1+5*x.^2)
f=sin(x)
df=diff(f)

nP=5 
xrange=[-pi,pi]
nnn=1600
X=[xrange(1):(xrange(2)-xrange(1))/(nP-1):xrange(2)] 
xxx=[xrange(1):(xrange(2)-xrange(1))/nnn:xrange(2)];

F=eval(subs(f,x,sym(X)))
DF=eval(subs(df,x,sym(X)))

n=length(X)

spalvos=['b','r','g','c','m','k','b','r','g','c','m','k'];
figure(1); hold on; grid on; title('Ermito daugianariai U');axis equal
figure(2); hold on; grid on; title('Ermito daugianariai V');axis equal
F
fff=0;
legU{1}='duoti interpoliavimo mazgai';legV{1}='duotos isvestiniu reiksmes';
legU{2}='interpoliavimo abscises';legV{2}='interpoliavimo abscises';
figure(1);plot(X,F,'ko','LineWidth',2,'MarkerSize',8);
          plot(X,zeros(size(X)),'k*');
          hU=legend(legU);
figure(2);plot(X,DF,'ko','LineWidth',2,'MarkerSize',8);
          plot(X,zeros(size(X)),'k*');hV=legend(legV);
          
for j=1:n
    [U,V]=Hermite(X,j,xxx);
    legU{j+2}=sprintf('U%2d',j);
    legV{j+2}=sprintf('V%2d',j);
    figure(1);plot(xxx,U,[spalvos(j),'-']);delete(hU);hU=legend(legU);
    figure(2);plot(xxx,V,[spalvos(j),'-']);delete(hV);hV=legend(legV);
    fff=fff+U*F(j)+V*DF(j);
    pause
end
 legU



figure(3), hold on, grid on, axis equal
plot(xxx,eval(subs(f,x,xxx)),'b-');
plot(xxx,fff,'r-','LineWidth',2.5);
legend({['duotoji funkcija ',char(f)],[' Ermito aproksimacija per ',sprintf('%d',n),' taskus']})
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
DL=0; %DL israiskos skaitiklis
for i=1:n % ciklas per atmetamus narius
    if i==j, continue, end 
    Lds=1;
    for k=1:n  
        if k ~= j && k ~= i , Lds=Lds.*(x-X(k)); end
    end
    DL=DL+Lds;
end
    Ldv=1;   %DL israiskos vardiklis 
    for k=1:n 
        if k ~= j, Ldv=Ldv.*(X(j)-X(k)); end 
    end
DL=DL/Ldv;

return
end
