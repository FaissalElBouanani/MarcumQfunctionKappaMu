function [R] = asymp(kapa,mu,m,sdb)
a=sqrt(2*(1-sqrt(0.5)));
b=sqrt(2*(1+sqrt(0.5)));
dirac=sqrt(b/a);
eta=(1-dirac^2)/dirac;
sigm=b-a;
sigma=b+a;
for j=1:length(sdb)
  snr=10^(sdb(j)/10); 
  lam=((m^m)*(mu*(1+kapa))^mu)/((mu*kapa+m)^m)/gamma(mu)/(snr^mu);
  if snr<1
      A=[0.1786 0.7564 dirac^2/(1-dirac^2)];
      B=[2.903 0.1307 0];
      C=[0.1786 0.7564 -0.5];
  else
      if snr>=0.1&snr<8
      A=[0.3798 0.6183 dirac^2/(1-dirac^2)];
      B=[1.895 7.93*10^(-4) 0];
      C=[0.3798 0.6183 -0.5];
      else
      A=[0.005206 0.6146 dirac^2/(1-dirac^2)];
      B=[0.2764 5.593*10^(-4) 0];
      C=[0.005206 0.6146 -0.5];
      end
  end
  omega=(mu^2)*kapa*(1+kapa)/(mu*kapa+m)/snr;
  
J1=A.*(1/2/sqrt(pi)).*((2./(2*B+sigm^2)).^(mu))*(gamma(mu+0.5)/(mu)).*hypergeom([0.5,mu],[mu+1],2*B./(sigm^2+2*B));
J2=A.*(1/2/sqrt(pi)).*((2./(2*B+sigma^2)).^(mu))*(gamma(mu+0.5)/(mu)).*hypergeom([0.5,mu],[mu+1],2*B./(sigma^2+2*B));
J3=C.*(gamma(mu)./((2+B).^mu)).*hypergeom([(mu)/2,(mu+1)/2],[1],2./(2+B).^2);
R(j)=lam*eta*(sum(J1)-sum(J2)+(1/eta)*sum(J3));
end
semilogy(sdb,R,'r');
end