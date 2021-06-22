% Parametrinis interpoliavimas_Ermito_splainais
% kiekvienas intervalas tarp tasku interpoliuojamas
% 2 eiles Ermito daugianariais, nustatant isvestiniu reiksmes pagal 
% skaitinio diferencijavimo formules
% valdomos tasku padetys ir isvestiniu reiksmes


function main
clc,close all

% spalva={'k-'};
spalva={'b-','r-','g-','m-','c-','k-'};

hL=[];
% X=[-1  0.5  1  1  -1 -1 ]
% Y=[ 1  0.5  1 -1  -1  1 ]

X=[0  0.2  3  4  9  10.5  14 14.5]
Y=[0.2 1 1.2 3  3  1.2  1  0.2]

nP=length(X);  % interpoliavimo mazgu skaicius
t(1)=0; for i=2:nP, t(i)=t(i-1)+norm([X(i) Y(i)]-[X(i-1) Y(i-1)]);  end  % parametro reiksmes
t

DX=Akima(t,X),DY=Akima(t,Y)
dist=(max(X)-min(X))/15 % posukiu rankeneles ilgis

f=figure(1);axis([min(X)-3*dist,max(X)+3*dist,min(Y)-3*dist,max(Y)+3*dist]);axis equal;hold on;

% vaizduojame duotus taskus
for i=1:nP, 
    h(i)=plot(X(i), Y(i),'ko','ButtonDownFcn',@startDragFcn,'MarkerSize',10);
    % kas atliekama paspaudus peles klavisa, nurodoma funkcijoje startDragFcn
    % tasku objektu valdikliai issaugomi masyve h
end

% vaizduojame liestines
for i=1:nP, 
    ang=atan2(DY(i),DX(i));
    ht(i)=plot(X(i)+dist*sin(ang), Y(i)-dist*cos(ang),'c*','ButtonDownFcn',@startDragFcn_t,'MarkerSize',10);
    % kas atliekama paspaudus peles klavisa, nurodoma funkcijoje startDragFcn_t
    % liestiniu objektu valdikliu taskai issaugomi masyve ht
end

set(f,'WindowButtonUpFcn',@stopDragFcn); % kas atliekama atleidus peles klavisa, nurodoma funkcijoje stopDragFcn

Ermito_splainu_parametrinis_interpoliavimas(X,DX,Y,DY,t);  % interpoliuojame pagal ivestus taskus ir 
                                                           % nubraizome pradine kreive  
                                                   
%************************************************************************   
% Toliau programa laukia pertraukimo nuo peles klaviso, kuris inicijuja 
% startDragFcn arba stopDragFcn vykdyma. Jos savo ruoztu peles judesi susieja arba atsieja 
% su draggingFcn
%************************************************************************
                                                   
%-----------  vidines funkcijos ------------------
%  jos aprasomos anksciau, nei sutinkamas pagrindines funkcijos "end",
%  todel visi pagrindineje funkcijoje naudojami kintamieji matomi taip pat
%  ir vidinese funkcijose

function startDragFcn(varargin)  % apraso, kas atliekama, kai paspaudziamas kairys peles klavisas   
    set(gcf, 'WindowButtonMotionFcn',@draggingFcn); % nurodo funkcija, kuria reikia nuolat kviesti pelei judant
end

function startDragFcn_t(varargin) % apraso, kas atliekama, kai paspaudziamas kairys peles klavisas
    set(gcf, 'WindowButtonMotionFcn',@draggingFcn_t); % nurodo funkcija, kuria reikia nuolat kviesti pelei judant
end

function draggingFcn(varargin)  % apraso, kas atliekama, kai pakinta pele valdomo objekto padetis
    iob=find(gco == h); % kelintas taskas judinamas
    x=get(h(iob),'xData');y=get(h(iob),'yData');  % tasko koordinates
    xt=get(ht(iob),'xData');yt=get(ht(iob),'yData'); % liestines valdiklio koordinates
    pt=get(gca,'Currentpoint');  % perskaitoma nauja padetis
    set(gco,'xData',pt(1,1),'yData',pt(1,2)); % pakeiciamos objekto koordinates
    X(iob)=pt(1,1);   % naujos tasko koordinates i masyva
    Y(iob)=pt(1,2);
    set(ht(iob),'xData',xt+pt(1,1)-x,'yData',yt+pt(1,2)-y); % pakeiciamos objekto koordinates
        % kvieciame savo sukurta funkcija interpoliuojanciai kreivei apskaiciuoti: 
    Ermito_splainu_parametrinis_interpoliavimas(X,DX,Y,DY,t);
end

function draggingFcn_t(varargin)    
    % apraso, kas atliekama, kai pakinta pele valdomo liestines valdiklio padetis
    iob=find(gco == ht); % kelinto tasko liestines valdiklis judinamas
    x=get(h(iob),'xData');y=get(h(iob),'yData');  % tasko koordinates
    xt=get(ht(iob),'xData');yt=get(ht(iob),'yData'); % liestines valdiklio koordinates
    ang0=atan2(yt-y,xt-x); % anksciau buves valdiklio rankeneles kampas
    pt=get(gca,'Currentpoint');  % perskaitoma nauja padetis
    dx=pt(1,1)-x; dy=pt(1,2)-y;  % valdiklio vektoriaus projekcijos 
    normd=norm([dx dy]); dx=dx/normd*dist;dy=dy/normd*dist;
    set(gco,'xData',x+dx,'yData',y+dy); % pakeiciamos objekto koordinates
    ang=atan2(dy,dx);  % naujas valdiklio rankeneles kampas
    DY(iob)=DY(iob)+sign(DX(iob))*tan(ang-ang0);
        % kvieciame savo sukurta funkcija interpoliuojanciai kreivei apskaiciuoti: 
    Ermito_splainu_parametrinis_interpoliavimas(X,DX,Y,DY,t);
end 

function stopDragFcn(varargin) % apraso, kas atliekama, kai atleidziamas kairys peles klavisas
    set(gcf, 'WindowButtonMotionFcn','');% nurodo, kad atleidus peles klavisa peles judejimas nebeturi kviesti funkcijos
end

function Ermito_splainu_parametrinis_interpoliavimas(X,DX,Y,DY,t);
    nP=length(X); % interpoliavimo tasku skaicius
    if ~isempty(hL), delete(hL); end

    for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
        nnn=100; 
        ttt=[t(iii):(t(iii+1)-t(iii))/nnn:t(iii+1)];
        fffX=0;fffY=0;
        for j=1:2
            [U,V]=Hermite(t(iii:iii+1),j,ttt);
            fffY=fffY+U*Y(iii+j-1)+V*DY(iii+j-1);
            fffX=fffX+U*X(iii+j-1)+V*DX(iii+j-1);
        end
        lsp=length(spalva);
        hL(iii)=plot(fffX,fffY,spalva{max(1,iii-floor(iii/lsp)*lsp)},'LineWidth',2.5);
    end %-----------------ciklas per intervalus pabaiga
return
end

function [U,V]=Hermite(X,j,x)   % Ermito daugianariai
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
        for k=1:n, if k ~= j && k ~= i , Lds=Lds.*(x-X(k)); end, end
        DL=DL+Lds;
    end
    Ldv=1;   %DL israskos vardiklis 
    for k=1:n, if k ~= j, Ldv=Ldv.*(X(j)-X(k)); end ,    end
    DL=DL/Ldv;
return
end

function DY=Akima(X,Y)
% Isvestiniu reiksmiu interpoliavimo taskuose nustatymas pagal skaitinio integravimo formules
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