function [h] = asympMPSK(kapa,mu,m,M,sdb)
g=(log(M)/log(2))*(sin(pi/M))^2;
A=[(7*M-8)/(48*M)  1/8  1/8  1/8  (M-2)/(12*M) (M-2)/(6*M) (M-2)/(6*M)];
dirac=[g 2*g  20*g/3  20*g/17  (g/(cos((M-2)*pi/(2*M))).^2)  (g/(cos((M-2)*pi/(6*M))).^2) (g/(cos((M-2)*pi/(3*M))).^2)];
%%%%%%%%%%%ù
lam=((m^m)*(mu*(1+kapa))^mu)/((mu*kapa+m)^m);
v=mu*(1+kapa);
omega=(mu^2)*kapa*(1+kapa)/(mu*kapa+m);
for k=1:length(sdb)
    snr=10^(sdb/10);
    lam=lam/(snr^mu);v=v/snr;omega=omega/snr;
    for jj=1:7
        
        s(jj)=(A(jj)*(dirac(jj))^(-mu));%
        
    end
h=lam*sum(s);
end

end
