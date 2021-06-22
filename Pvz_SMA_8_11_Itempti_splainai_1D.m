% Interpoliavimas itemptais splainais
 
function Splainu_interpoliavimas
clc,close all
figure(1), hold on, grid on, axis equal
syms  f x 

% f=sin(x)  % duotoji funkcija
f=1./(1+5*x.^2)
df=diff(f)
nP=5 % interpoliavimo tasku skaicius
sgm=[-1 -5 -5 -1]*1
% itempimai segmentuose

xrange=[-pi,pi]
X=[xrange(1):(xrange(2)-xrange(1))/(nP-1):xrange(2)] 
Y=eval(subs(f,sym(x),sym(X)))
plot(X,Y,'ko');

DDF=itempto_splaino_koeficientai(X,Y,sgm,0)
for iii=1:nP-1  %------  ciklas per intervalus tarp gretimu tasku
    nnn=100;
    [S,sss]=itemptas_splainas(X(iii:iii+1),Y(iii:iii+1),DDF(iii:iii+1),sgm(iii),nnn);
    plot(sss,eval(subs(f,sym(x),sym(sss))),'b-');
    plot(sss,S,'k-','LineWidth',2,'MarkerSize',8)
end %-----------------ciklo per intervalus pabaiga
legend({['duotoji funkcija',char(f)],[sprintf('Itempti splainai %d intervaluose   sgm=',nP-1),sprintf('%g ',sgm)]});

return
end


function DDF=itempto_splaino_koeficientai(X,Y,sg,iopt)
% apskaiciuojamos antros isvestines splaino mazguose
% iopt=1 - periodinis splainas
% ipot=0 - splaino galuose "sarnyrai" 
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


function [S,sss]=itemptas_splainas(X,Y,DDF,sgm,nnn)
% splaino intervale tarp dvieju tasku apskaiciavimas
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
