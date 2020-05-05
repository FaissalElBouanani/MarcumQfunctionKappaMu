function [ ] = diverst
x=[1:0.5:60];kapa=5;M=31;
m=4.7;
%%%%%%%%%%%%%%%%%
  mu=1.9 ;
 for j=1:length(x)
    sdb=x(j);
R=anal2(kapa,mu,m,sdb,M);
y1(j)=abs(log(R)./log(10.^(sdb./10)))
end
semilogy(x,y1,'k');%
hold on
%%%%%%%%%%%%%%%%%%%
 for j=1:length(x)
    sdb=x(j);
R=anal2(kapa,mu,m,sdb,M);
y1(j)=abs(log(R)./log(10.^(sdb./10)))
end
semilogy(x,y1,'r');%
hold on
