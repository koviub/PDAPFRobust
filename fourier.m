function [vs]=fourier(list)

n=length(list);
ur=zeros(n,n);
vs=zeros(1,n);
for ii=1:n
    
    rr=1:n;
    ur(ii,:)=list.*exp(2*pi*1i/n*(ii-1).*(rr-1));
    vs(ii)=1/sqrt(n).*sum(ur(ii,:));
    
end
end