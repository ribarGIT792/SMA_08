
% Itemptu_splainu_interpoliavimas_parametrinis
% Pele valdomos interpoliavimo tasku padetys

function main
clc,close all,clear all 
hL=[];  % splaino intervalu objekto rodykles
segment=[];  % relaksuojamo intervalo numeris

f=figure; hold on; grid on
iopt=0;
X=[-1  0.5  1  1  -1 -1 ]
Y=[ 1  0.5 1 -1  -1  1 ]
% X=[0  0.2  3  4  9  10.5  14 14.5]
% Y=[0.2 1 1.2 3  3  1.2  1  0.2]
nP=length(X);
t(1)=0; for i=2:nP, t(i)=t(i-1)+norm([X(i) Y(i)]-[X(i-1) Y(i-1)]);  end
t

sg(1:nP-1)=1  % itempimai
dd=(max(X)-min(X))/2;
figure(1);axis([min(X)-dd,max(X)+dd,min(Y)-dd,max(Y)+dd]);axis equal;hold on;

% vaizduojame duotus taskus
for i=1:nP,  h(i)=plot(X(i), Y(i),'ko','ButtonDownFcn',@startDragFcn,'MarkerSize',10); end 
    % kas atliekama paspaudus peles klavisa, nurodoma funkcijoje startDragFcn
    % tasku objektu valdikliai issaugomi masyve h

set(f,'WindowButtonUpFcn',@stopDragFcn); % kas atliekama atleidus peles klavisa, nurodoma funkcijoje stopDragFcn

itemptu_splainu_parametrinis_interpoliavimas(X,Y,t,sg);  % interpoliuojame pagal ivestus taskus ir 
                                                   % nubraizome pradine kreive  
                                                   
%************************************************************************   
% Toliau programa laukia pertraukimo nuo peles klaviso, kuris inicijuoja 
% startDragFcn arba stopDragFcn vykdyma. Jos savo ruoztu peles judesi susieja arba atsieja 
% su draggingFcn
%************************************************************************
                                                   
%-----------  vidines funkcijos ------------------
%  jos aprasomos anksciau, nei sutinkamas pagrindines funkcijos "end",
%  todel visi pagrindineje funkcijoje naudojami kintamieji matomi taip pat
%  ir vidinese funkcijose

function startDragFcn(varargin)   % apraso, kas atliekama, kai paspaudziamas kairys peles klavisas ties mazgu
    set(gcf, 'WindowButtonMotionFcn',@draggingFcn); % nurodo funkcija, kuria reikia nuolat kviesti pelei judant
end

function startTensFcn(varargin)
    % apraso, kas atliekama, kai paspaudziamas kairys peles klavisas ties
    % splaino segmentu: didinamas splaino segmento itempimas
    segment=find(gco == hL);
    sg(segment)=sg(segment)+1;
    itemptu_splainu_parametrinis_interpoliavimas(X,Y,t,sg);
    set(gcf, 'WindowButtonMotionFcn',@relaxingFcn); % nurodo funkcija, kuria reikia nuolat kviesti pelei judant
end

function relaxingFcn(varargin) % apraso, kas atliekama, kai pele juda; judinant pele segmentas relaksuojamas
    sg(segment)=max(sg(segment)-1,1);
    itemptu_splainu_parametrinis_interpoliavimas(X,Y,t,sg);
end

function draggingFcn(varargin)  % apraso, kas atliekama, kai pele juda
    pt=get(gca,'Currentpoint'); % perskaitoma nauja padetis
    set(gco,'xData',pt(1,1),'yData',pt(1,2)); % pakeiciamos objekto koordinates
    X(find(gco == h))=pt(1,1); Y(find(gco == h))=pt(1,2);
    % kvieciame savo sukurta funkcija interpoliuojanciai kreivei apskaiciuoti: 
    itemptu_splainu_parametrinis_interpoliavimas(X,Y,t,sg);
end

function stopDragFcn(varargin) % apraso, kas atliekama, kai atleidziamas kairys peles klavisas
    set(gcf, 'WindowButtonMotionFcn','');% nurodo, kad atleidus peles klavisa peles judejimas nebeturi kviesti funkcijos
end

function itemptu_splainu_parametrinis_interpoliavimas(X,Y,t,sg) 
    nP=length(X) % interpoliavimo tasku skaicius
    if ~isempty(hL), delete(hL); end
%     iopt=1;
    DDFX=itempto_splaino_koeficientai(t,X,sg,iopt);
    DDFY=itempto_splaino_koeficientai(t,Y,sg,iopt);
    
    fclose all,  % splaino vaizdavimo taskai irasomi i failus
    delete('carx.txt') , delete('cary.txt') ,
    fhx=fopen('carx.txt','a+'); fhy=fopen('cary.txt','a+');
    
    for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
        nnn=100;
        [SX,sss]=itemptas_splainas(t(iii:iii+1),X(iii:iii+1),DDFX(iii:iii+1),sg(iii),nnn);
        [SY,sss]=itemptas_splainas(t(iii:iii+1),Y(iii:iii+1),DDFY(iii:iii+1),sg(iii),nnn);
        hL(iii)=plot(SX,SY,'k-','LineWidth',2,'MarkerSize',8,'ButtonDownFcn',@startTensFcn);
        nn1=length(SX); if iii < nP-1, nn1=nn1-1; end
        fprintf(fhx,'%g ',SX(1:nn1)); 
        fprintf(fhy,'%g ',SY(1:nn1)); 
    end %-----------------ciklas per intervalus pabaiga
    fclose all
return
end

function DDF=itempto_splaino_koeficientai(X,Y,sg,iopt) % apskaiciuojamos antros isvestines splaino mazguose
% iopt=1 - periodinis splainas
    n=length(X);
    A=zeros(n);b=zeros(n,1);
    d=X(2:n)-X(1:(n-1));
    for i=1:n-2
        sg1=sg(i);sg2=sg(i+1);
        A(i,i:i+2)=[1/(sg1^2*d(i))-1/(sg1*sinh(sg1*d(i))), ...
                 cosh(sg1*d(i))/(sg1*sinh(sg1*d(i)))+cosh(sg2*d(i+1))/(sg2*sinh(sg2*d(i+1)))-1/(sg1^2*d(i))-1/(sg2^2*d(i+1)),...
                 1/(sg2^2*d(i+1))-1/(sg2*sinh(sg2*d(i+1)))];
        b(i)=(Y(i+2)-Y(i+1))/d(i+1)-(Y(i+1)-Y(i))/d(i);
    end
    if iopt == 0,   A(n-1,1)=1;A(n,n)=1;b(n-1)=0;b(n)=0;
    else,sg1=sg(1);sgd1=sg(1)*d(1);sgn1=sg(n-1);sgdn1=sg(n-1)*d(n-1);
         A(n-1,[1,2,n-1,n])=[-cosh(sgd1)/sinh(sgd1)/sg1+1/sgd1/sg1,1/sinh(sgd1)/sg1-1/sgd1/sg1, ...
                             1/sinh(sgdn1)/sgn1-1/sgdn1/sgn1, -cosh(sgdn1)/sinh(sgdn1)/sgn1+1/sgdn1/sgn1]; 
         A(n,[1,n])=[1,-1];  
         b(n-1)=(Y(1)-Y(2))/d(1)-(Y(n-1)-Y(n))/d(n-1);
    end
    DDF=A\b;
return
end


function [S,sss]=itemptas_splainas(X,Y,DDF,sgm,nnn) % splaino intervale tarp dvieju tasku apskaiciavimas
% nnn - vaizdavimo tasku skaicius
% S - splaino reiksmes
% sss - vaizdavimo abscises
    d=X(2)-X(1);
    sss=X(1):(X(2)-X(1))/(nnn-1):X(2);
    S=DDF(1)/sgm^2*sinh(sgm*(d-(sss-X(1))))/sinh(sgm*d)+...
        (Y(1)-DDF(1)/sgm^2)*(d-(sss-X(1)))/d + ...
        DDF(2)/sgm^2*sinh(sgm*((sss-X(1))))/sinh(sgm*d)+...
        (Y(2)-DDF(2)/sgm^2)*((sss-X(1)))/d;
return
end

end   % Sis end uzbaigia pagrindine funkcija