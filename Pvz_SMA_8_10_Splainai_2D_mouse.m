% Splainu_interpoliavimas_2D_parametrinis
% Pele valdomos interpoliavimo tasku padetys

function main
clc,close all,clear all 

spalva={'k-'};
% spalva={'b-','r-','g-','m-','c-','k-'};

hL=[]; % busimu objektu valdikliu masyvas
f=figure; hold on; grid on

X=[-1  0.5  1  1  -1 -1 ]
Y=[ 1  0.5  1 -1  -1  1 ]
% X=[-1  0  1  ]
% Y=[ 1  1  1  ]

global iper
iper=1;  % 0 - splainas laisvais galais 
         % 1 - periodinis splainas 

nP=length(X);
t(1)=0; for i=2:nP, t(i)=t(i-1)+norm([X(i) Y(i)]-[X(i-1) Y(i-1)]);  end
t

figure(1);axis([-3,3,-3,3]);axis equal;hold on;

% vaizduojame duotus taskus
for i=1:nP, 
    h(i)=plot(X(i), Y(i),'ko','ButtonDownFcn',@startDragFcn,'MarkerSize',10,'MarkerFaceColor','k');
    % kas atliekama paspaudus peles klavisa, nurodoma funkcijoje startDragFcn
    % tasku objektu valdikliai issaugomi masyve h
end

set(f,'WindowButtonUpFcn',@stopDragFcn); % kas atliekama atleidus peles klavisa, nurodoma funkcijoje stopDragFcn

splainu_parametrinis_interpoliavimas(X,Y,t);  % interpoliuojame pagal ivestus taskus ir 
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

function startDragFcn(varargin)
    % apraso, kas atliekama, kai paspaudziamas kairys peles klavisas
    set(gcf, 'WindowButtonMotionFcn',@draggingFcn); % nurodo funkcija, kuria reikia nuolat kviesti pelei judant
end

function draggingFcn(varargin)    
    % apraso, kas atliekama, kai pakinta pele valdomo objekto padetis
    
    pt=get(gca,'Currentpoint');  % perskaitoma nauja padetis
    set(gco,'xData',pt(1,1),'yData',pt(1,2)); % pakeiciamos objekto koordinates
    X(find(gco == h))=pt(1,1);   
    Y(find(gco == h))=pt(1,2);
    % kvieciame savo sukurta funkcija interpoliuojanciai kreivei apskaiciuoti: 
    splainu_parametrinis_interpoliavimas(X,Y,t);
    
end

function stopDragFcn(varargin)
    % apraso, kas atliekama, kai atleidziamas kairys peles klavisas
    set(gcf, 'WindowButtonMotionFcn','');% nurodo, kad atleidus peles klavisa peles judejimas nebeturi kviesti funkcijos
    
end

function splainu_parametrinis_interpoliavimas(X,Y,t)
nP=length(X) % interpoliavimo tasku skaicius
if ~isempty(hL), delete(hL); end
% 
% iper=1;  % 0 - splainas laisvais galais 
%          % 1 - periodinis splainas  
DDFX=splaino_koeficientai(t,X,iper);
DDFY=splaino_koeficientai(t,Y,iper);

for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku

nnn=100;
[SX,sss]=splainas(t(iii:iii+1),X(iii:iii+1),DDFX(iii:iii+1),nnn);
[SY,sss]=splainas(t(iii:iii+1),Y(iii:iii+1),DDFY(iii:iii+1),nnn);
lsp=length(spalva);
spalva{max(1,iii-floor(iii/lsp)*lsp)}
hL(iii)=plot(SX,SY,spalva{max(1,iii-floor(iii/lsp)*lsp)},'LineWidth',2,'MarkerSize',8)

end %-----------------ciklas per intervalus pabaiga
% splaino intervalu objektu valdikliai issaugomi masyve hL

return
end


function DDF=splaino_koeficientai(X,Y,iper)
% apskaiciuojamos antros isvestines splaino mazguose
% iopt=1 - periodinis splainas

n=length(X);
A=zeros(n);b=zeros(n,1);
d=X(2:n)-X(1:(n-1));
 for i=1:n-2
     A(i,i:i+2)=[d(i)/6, (d(i)+d(i+1))/3,d(i+1)/6];
     b(i)=(Y(i+2)-Y(i+1))/d(i+1)-(Y(i+1)-Y(i))/d(i);
 end
 
if iper == 0,  A(n-1,1)=1;A(n,n)=1;
else, A(n-1,[1,2,n-1,n])=[d(1)/3, d(1)/6, d(n-1)/6,d(n-1)/3]; 
      A(n,[1,n])=[1,-1];  
      b(n-1)=(Y(2)-Y(1))/d(1)-(Y(n)-Y(n-1))/d(n-1);
end

DDF=A\b;
 
return
end


function [S,sss]=splainas(X,Y,DDF,nnn)
% splaino intervale tarp dvieju tasku apskaiciavimas
% nnn - vaizdavimo tzku skaicius
% S - splaino reiksmes
% sss - vaizdavimo abscises
d=X(2)-X(1);
sss=X(1):(X(2)-X(1))/(nnn-1):X(2);
S=DDF(1)/2*(sss-X(1)).^2+(DDF(2)-DDF(1))/(6*d)*(sss-X(1)).^3+(sss-X(1))*((Y(2)-Y(1))/d-DDF(1)*d/3-DDF(2)*d/6) +Y(1);

return
end


end   % Sis end uzbaigia pagrindine funkcija
