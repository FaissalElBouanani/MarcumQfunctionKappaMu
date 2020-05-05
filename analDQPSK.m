function [RR] = analDQPSK(kapa,mu,m,sdb,M)
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
  xi=mu*(1+kapa)/snr+B;
  for k=0:M  
g(k+1)=sum(((2./(2*xi+sigm^2)).^mu).*A.*((pochhammer(m,k)/pochhammer(mu,k))/(factorial(k))).*((2*omega./(2*xi+sigm^2)).^(k))*(gamma(mu+k+0.5)/(mu+k)).*hypergeom([0.5,mu+k],[mu+k+1],2*xi./(sigm^2+2*xi)));
end
I1(j)=(1/2/sqrt(pi))*sum(g);
   for k=0:M  
g(k+1)=sum(((2./(2*xi+sigma^2)).^mu).*A.*((pochhammer(m,k)/pochhammer(mu,k))/(factorial(k))).*((2*omega./(2*xi+sigma^2)).^(k))*(gamma(mu+k+0.5)/(mu+k)).*hypergeom([0.5,mu+k],[mu+k+1],2*xi./(sigma^2+2*xi)));
end
I2(j)=(1/2/sqrt(pi))*sum(g);
for k=0:M
g(k+1)=sum(C.*((1./(2+xi)).^mu).*(pochhammer(m,k)./pochhammer(mu,k))./(factorial(k)).*((omega./(2+xi)).^k).*(gamma(mu+k)).*hypergeom([(mu+k)/2,(mu+k+1)/2],[1],2./(2+xi).^2));
end
I3(j)=sum(g);
R(j)=lam*eta*(I1-I2+(1/eta)*I3);
end


end