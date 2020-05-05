function [rr] = exactMPSK(kapa,mu,m,M,sdb)
lam=((m^m)*(mu*(1+kapa))^mu)/((mu*kapa+m)^m)/gamma(mu);
v=mu*(1+kapa);
omega=(mu^2)*kapa*(1+kapa)/(mu*kapa+m);
A=(log(M)/log(2))*(sin(pi/M))^2;
t=[0.01:0.1:50];
for j=1:length(t)
    x=t(j);
    teta=[0:0.01:(M-1)*pi/M];
    g=exp(-A*x./(sin(teta)).^2);
    pexact(j)=(1/pi)*trapz(teta,g);
end
z=t;
for k=1:length(sdb)
    snr=10^(sdb(k)/10);
    lam=lam/(snr^mu);v=v/snr;omega=omega/snr;
    f=lam.*z.^(mu-1).*exp(-v*z).*hypergeom([m],[mu],omega*z);%
    ff=f.*pexact;
    rr(k)=trapz(z,ff);
end

end