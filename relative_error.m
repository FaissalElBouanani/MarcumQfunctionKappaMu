function [] = relative_error
x=[0.1:0.1:16];

for j=1:length(x)
    %%%%%%%Calcule exact
    a=sqrt(2*x(j)*(1-sqrt(0.5)));
    b=sqrt(2*x(j)*(1+sqrt(0.5)));%2*sqrt(x(j));%
    %a=ss*b;
    pexact(j)=marcumq(a,b)-0.5*besseli(0,a*b)*exp(-0.5*(a^2+b^2));
       %barzic Approximations BER5, BER6 and BER7
     L1(j)= besseli(0,a*b)*(sqrt(pi/2)*(b/exp(a*b))*erfc((b-a)/sqrt(2))-0.5*exp(-0.5*(a^2+b^2)));
     L2(j)=besseli(0,a*b)*((sqrt(0.5*pi)*b/(exp(a*b)-exp(-a*b)))*(erfc((b-a)/sqrt(2))-erfc((b+a)/sqrt(2)))-0.5*exp(-0.5*(a^2+b^2)));
     U1(j)= besseli(0,a*b)*(sqrt(pi/2)*(a/exp(a*b))*erfc((b-a)/sqrt(2))+0.5*exp(-0.5*(a^2+b^2)));
     U2(j)=besseli(0,a*b)*((sqrt(0.5*pi)*a/(exp(a*b)+exp(-a*b)))*(erfc((b-a)/sqrt(2))-erfc((b+a)/sqrt(2)))+0.5*exp(-0.5*(a^2+b^2)));
     U3(j)=besseli(0,a*b)*((sqrt(0.5*pi)*a/(exp(a*b)+3))*erfc((b-a)/sqrt(2))+0.5*exp(-0.5*(a^2+b^2)));
     %%%%%%%%%%%BER5
     if x(j)<=1
         a5=0.65*(x(j))^(0.25);
     else
         a5=0.5+0.5*sqrt(2)*1.1*exp(-0.5*pi/sqrt(x(j)))/(x(j)^(1.5));
     end
     ber5(j)=a5*L1(j)+(1-a5)*U1(j);
     
     %%%%%%%%%%BER6
     if x(j)<1
        a6=0.25*exp(-(x(j)^2)/2.9)+0.5;
    else
    if x(j)<5&x(j)>=1
       a6=((exp(-1/(2*x(j)+1)))/((x(j)+0.5)^(1.5)))*1.15*sqrt(0.5/pi)+0.5;
    else
        a6=0.5+(0.65/pi)/(1+x(j));
        end
     end
    ber6(j)=a6*L2(j)+(1-a6)*U2(j);
     %%%%%%BER7
     if x(j)<1
        a7=0.95*(1-x(j))^2;
    else
    if x(j)<=8&x(j)>=1
       a7=0.5-1.4*(exp(-x(j)^1.2)+0.02);% for x entre 0.1 et 1
    else
        a7=0.5+1/(5.2*x(j));
        end
    end
    ber7(j)=a7*L2(j)+(1-a7)*U3(j);
    %Poposed approximation
     a1=sqrt(2*(1-sqrt(0.5)));
    b1=sqrt(2*(1+sqrt(0.5)));
    UPP(j)=sqrt(a1/b1)*(-qfunc((b1+a1)*sqrt(x(j)))+qfunc((b1-a1)*sqrt(x(j))))+0.5*besseli(0,sqrt(2)*x(j)).*exp(-2*x(j));%0.5*(1/sqrt(2*pi*sqrt(2)*x(j)))*exp(sqrt(2)*x(j)).*exp(-2*x(j));%0.5*besseli(0,sqrt(2)*x(j)).*exp(-2*x(j));%
    LR(j)=sqrt(b1/a1)*(-qfunc((b1+a1)*sqrt(x(j)))+qfunc((b1-a1)*sqrt(x(j))))-0.5*besseli(0,sqrt(2)*x(j)).*exp(-2*x(j));%0.5*(1/sqrt(2*pi*sqrt(2)*x(j)))*exp(sqrt(2)*x(j)).*exp(-2*x(j));%0.5*besseli(0,sqrt(2)*x(j)).*exp(-2*x(j));
    %%%%%%%%%%%%%%%%%%%%%%our Proposed
if x(j)<1
        aa=0.1786*exp(-2.903*x(j))+0.7564*exp(-0.1307*x(j));
else
    if x(j)>=1&x(j)<=8
   aa=0.3798*exp(-1.895*x(j))+0.6183*exp(-0.000793*x(j)); 
    else
     aa=0.005206*exp(-0.2764*x(j))+0.6146*exp(-5.593*10^(-5)*x(j));
    end    
    end
Tupp(j)=aa*UPP(j)+(1-aa)*LR(j);
end
 er1=abs(Tupp-pexact)./pexact;
 er2=abs(ber5-pexact)./pexact;
 er3=abs(ber6-pexact)./pexact;
 er4=abs(ber7-pexact)./pexact;
   %ber6
%  ber7
% % %  
  %semilogy(sdb,er1,'r',sdb,er3,'b',sdb,er2,'k',sdb,er4,'g')
  %pexact
%   Tupp
%   ber6
%   ber7
%  hold on
  %plot(sdb,Tupp,'r',sdb,ber6,'b',sdb,ber7,'k',sdb,pexact,'k--')
  
%   hold on
%   semilogy(sdb,pexact,'k',sdb,barzic,'b',sdb,TUPP,'r');%,sdb,U2,'k--',sdb,U3,'b--');%,sdb,TUPP,'r')
% hold on

%end
end
