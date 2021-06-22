% Interpoliavimas splainais
% 
function Splainu_interpoliacija
clc,close all

figure(1), hold on, grid on, axis equal
spalva={'k-'};
% spalva={'b-','r-','g-','m-','c-','k-'};
syms  f x 

f=sin(x)  % duotoji funkcija
f=1./(1+5*x.^2)
df=diff(f)
nP=15 % interpoliavimo tasku skaicius
xrange=[-pi,pi]
X=[xrange(1):(xrange(2)-xrange(1))/(nP-1):xrange(2)] 
Y=eval(subs(f,sym(x),sym(X)));
plot(X,Y,'ko');

DDF=splaino_koeficientai(X,Y,1)


for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku

nnn=100;
[S,sss]=splainas(X(iii:iii+1),Y(iii:iii+1),DDF(iii:iii+1),nnn);
plot(sss,eval(subs(f,sym(x),sym(sss))),'b-');
lsp=length(spalva);
plot(sss,S,spalva{max(1,iii-floor(iii/lsp)*lsp)},'LineWidth',2,'MarkerSize',8)
end %-----------------ciklas per intervalus pabaiga
legend({['duotoji funkcija ',char(f)],sprintf('Splainai %d intervaluose',nP-1)});
return
end




function DDF=splaino_koeficientai(X,Y,iopt)
% apskaiciuojamos antros isvestines splaino mazguose
% iopt=1 - periodinis splainas

n=length(X);
A=zeros(n);b=zeros(n,1);
d=X(2:n)-X(1:(n-1))
 for i=1:n-2
     A(i,i:i+2)=[d(i)/6, (d(i)+d(i+1))/3,d(i+1)/6];
     b(i)=(Y(i+2)-Y(i+1))/d(i+1)-(Y(i+1)-Y(i))/d(i);
 end
 
if iopt == 0,  A(n-1,1)=1;A(n,n)=1;
else, A(n-1,[1,2,n-1,n])=[d(1)/3, d(1)/6, d(n-1)/6,d(n-1)/3]; 
      A(n,[1,n])=[1,-1];  
      b(n-1)=(Y(2)-Y(1))/d(1)-(Y(n)-Y(n-1))/d(n-1);
end

DDF=A\b;
 
return
end


function [S,sss]=splainas(X,Y,DDF,nnn)
% splaino intervale tarp dvieju tasku apskaiciavimas
% nnn - vaizdavio tzku skaicius
% S - splaino reiksmes
% sss - vaizdavimo abscises
d=X(2)-X(1);
sss=X(1):(X(2)-X(1))/(nnn-1):X(2);
S=DDF(1)/2*(sss-X(1)).^2+(DDF(2)-DDF(1))/(6*d)*(sss-X(1)).^3+(sss-X(1))*((Y(2)-Y(1))/d-DDF(1)*d/3-DDF(2)*d/6) +Y(1);

return
end
