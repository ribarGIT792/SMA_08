% Interpoliavimas_Ermito_splainais
% kiekvienas intervalas tarp tasku interpoliuojamas
% 2 eiles Ermito daugianariais, nustatant isvestiniu reiksmes pagal 
% skaitinio diferencijavimo formules

function main
clc,close all
hL=[];

X=[0 1 2 3 4 5 6 7 8]
Y=sin(1./(1+2*(X-4).^2))
DY=Akima(X,Y)
nP=length(X) % interpoliavimo tasku skaicius

f=figure(1),hold on, axis equal, grid on
% vaizduojame duotus taskus
for i=1:nP, 
    h(i)=plot(X(i), Y(i),'ko','ButtonDownFcn',@startDragFcn,'MarkerSize',10);
    % kas atliekama paspaudus peles klavisa, nurodoma funkcijoje startDragFcn
    % tasku objektu valdikliai issaugomi masyve h
end

set(f,'WindowButtonUpFcn',@stopDragFcn); % kas atliekama atleidus peles klavisa, nurodoma funkcijoje stopDragFcn

Ermito_splainu_interpoliavimas(X,Y,DY);  % interpoliuojame pagal ivestus taskus ir nubraizome pradine kreive  
                                                   
%************************************************************************   
% Toliau programa laukia pertraukimo nuo peles klaviso, kuris inicijuja 
% startDragFcn arba stopDragFcn vykdyma. Jos savo ruoztu peles judesi susieja arba atsieja 
% su draggingFcn
%************************************************************************

                                                   
%-----------  vidines funkcijos ------------------
%  jos aprasomos anksciau, nei sutinkamas pagrindines funkcijos "end",
%  todel visi pagrindineje funkcijoje naudojami kintamieji matomi taip pat
%  ir vidinese funkcijose

function startDragFcn(varargin) % apraso, kas atliekama, kai paspaudziamas kairys peles klavisas
     set(gcf, 'WindowButtonMotionFcn',@draggingFcn); % nurodo funkcija, kuria reikia nuolat kviesti pelei judant
end

function draggingFcn(varargin)  % apraso, kas atliekama, kai pele juda
    pt=get(gca,'Currentpoint');  % perskaitoma nauja padetis
    set(gco,'xData',pt(1,1),'yData',pt(1,2)); % pakeiciamos objekto koordinates
    X(find(gco == h))=pt(1,1); Y(find(gco == h))=pt(1,2);
    % kvieciame savo sukurta funkcija interpoliuojanciai kreivei apskaiciuoti: 
    Ermito_splainu_interpoliavimas(X,Y,DY);
 end

function stopDragFcn(varargin) % apraso, kas atliekama, kai atleidziamas kairys peles klavisas
    set(gcf, 'WindowButtonMotionFcn','');% nurodo, kad atleidus peles klavisa peles judejimas nebeturi kviesti funkcijos
end

function Ermito_splainu_interpoliavimas(X,Y,DY);
    nP=length(X) % interpoliavimo tasku skaicius
    if ~isempty(hL), delete(hL); end
    for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
        nnn=100;
        xxx=[X(iii):(X(iii+1)-X(iii))/nnn:X(iii+1)];
        fff=0;
        for j=1:2
            [U,V]=Hermite(X(iii:iii+1),j,xxx);
            fff=fff+U*Y(iii+j-1)+V*DY(iii+j-1);
        end
        hL(iii)=plot(xxx,fff,'r-','LineWidth',2.5);
    end %-----------------ciklo per intervalus pabaiga
return
end

function [U,V]=Hermite(X,j,x)  % Ermito daugianariai
    L=Lagrange(X,j,x); DL=D_Lagrange(X,j,X(j));
    U=(1-2*DL.*(x-X(j))).*L.^2;
    V=(x-X(j)).*L.^2;
return
end

function L=Lagrange(X,j,x)  % Lagranzo daugianaris
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
        for k=1:n, if k ~= j && k ~= i , Lds=Lds.*(x-X(k)); end,  end
        DL=DL+Lds;
       end
    Ldv=1;   %DL israskos vardiklis 
    for k=1:n 
        if k ~= j, Ldv=Ldv.*(X(j)-X(k)); end 
    end
    DL=DL/Ldv;
return
end

function DY=Akima(X,Y) % Isvestiniu reiksmiu interpoliavimo taskuose nustatymas pagal skaitinio integravimo formules
    n=length(X);
    fnk=inline('(2*x-xi-xip1)/((xim1-xi)*(xim1-xip1))*yim1+(2*x-xim1-xip1)/((xi-xim1)*(xi-xip1))*yi+(2*x-xim1-xi)/((xip1-xim1)*(xip1-xi))*yip1')
    for i=1:n
        if i == 1,xim1=X(1);xi=X(2);xip1=X(3); yim1=Y(1);yi=Y(2);yip1=Y(3);DY(i)=fnk(xim1,xi,xim1,xip1,yi,yim1,yip1);
        elseif i == n, xim1=X(n-2);xi=X(n-1);xip1=X(n); yim1=Y(n-2);yi=Y(n-1);yip1=Y(n); DY(n)=fnk(xip1,xi,xim1,xip1,yi,yim1,yip1);
        else, xim1=X(i-1);xi=X(i);xip1=X(i+1); yim1=Y(i-1);yi=Y(i);yip1=Y(i+1); DY(i)=fnk(xi,xi,xim1,xip1,yi,yim1,yip1);
        end
    end
return
end


end   % Sis end uzbaigia pagrindine funkcija
