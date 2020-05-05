function [rr] = exact(kapa,mu,m,sdb)
a=sqrt(2*(1-sqrt(0.5)));
b=sqrt(2*(1+sqrt(0.5)));
dirac=sqrt(b/a);
eta=(1-dirac^2)/dirac;
lam=((m^m)*(mu*(1+kapa))^mu)/((mu*kapa+m)^m)/gamma(mu);
v=mu*(1+kapa);
omega=(mu^2)*kapa*(1+kapa)/(mu*kapa+m);
t=rand(1,5000);
       x=-log(t);
      pexact=marcumq(a*sqrt(x),b*sqrt(x))-0.5*besseli(0,sqrt(2)*x).*exp(-2*x);
for k=1:length(sdb)
        snr=10^(sdb(k)/10);
    lam=lam/(snr^mu);v=v/snr;omega=omega/snr;
    f=(1./t).*lam.*x.^(mu-1).*exp(-v*x).*hypergeom([m],[mu],omega*x);%
     ff=f.*pexact;
     rr(k)=mean(ff)  
end
end