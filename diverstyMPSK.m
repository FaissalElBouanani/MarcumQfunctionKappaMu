function [ ] = diverstMPSK

kapa=5;M=8;m=4.7;
%%%%%%%%%%%%%%%
  mu=1.9 
 for j=1:length(x)
    sdb=x(j);
h=analMPSK(kapa,mu,m,M,sdb)
y1(j)=abs(log(h)./log(10.^(sdb./10)));
end
semilogy(x,y1,'k--');%
hold on
%%%%%%%%%%%%%%%%%
 mu=2.9;
 for j=1:length(x)
    sdb=x(j);
h=analMPSK(kapa,mu,m,M,sdb)
y1(j)=abs(log(h)./log(10.^(sdb./10)));
end
semilogy(x,y1,'k--');%
hold on