function res = Robustness(input,full)
syms var1 var2 real
syms s

var=[var1 var2];

param=str2num(get(input.edit6,'String'));
vari=str2num(get(input.edit7,'String'));
weights=str2num(get(input.edit16,'String'));

x0=str2double(get(input.edit8,'String'));                          % origin on the complex plane
Rad=str2double(get(input.edit10,'String'));  

% structured perturbation on var1 var2 
F=str2func(get(input.edit1,'String'));
covar=coeffs(F(s,var,param),var);
Dc=@(x)double(F(x,vari,param));
% w1=@(x)double(subs([weights(1)*vari(1)*covar(3);weights(2)*vari(2)*covar(2);weights(3)*param(4)*s^2*exp(-param(2)*s)],'s',x));
% w1=@(x)double([weights(1)*vari(1);weights(2)*vari(2)*x;weights(3)*param(4)*x^2].*exp(-param(2)*x));% pda weight function
w1=@(x)double([weights(1)*vari(1)*(1 - (exp(-(x*(param(2) + param(2)*(param(3) + 1)))/2)*sinh((x*(param(2) - param(2)*(param(3) + 1)))/2)*(2*x*sinh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1)) + 2*(3/2*9.81/param(1))^(1/2)*cosh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1))))/(3/2*9.81/param(1))^(1/2));...
    weights(2)*vari(2)*(x - (exp(-(x*(param(2) + param(2)*(param(3) + 1)))/2)*sinh((x*(param(2) - param(2)*(param(3) + 1)))/2)*(2*(3/2*9.81/param(1))*sinh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1)) + 2*(3/2*9.81/param(1))^(1/2)*x*cosh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1))))/(3/2*9.81/param(1))^(1/2));...
    weights(3)*param(4)*x^2.*exp(-param(2)*x)]);% pf weight function

W=@(x)w1(x)./Dc(x);

re=-Rad:Rad/50:Rad/2;
im=-Rad:Rad/50:Rad;
[R,I]=meshgrid(re,im);
lambda=R+1i*I+x0;
if full
    % calculate the pseudo spectrum over complex domain
    r=SpectralRadius(W,lambda,'real');
else
    r=zeros(size(R));
end
% caulate the stability radius
lamb=1i*(-0:Rad/100:Rad);
[~,Ind]=min(SpectralRadius(W,lamb,'real'));
            if Ind==1
                Ind=2;
            elseif Ind==length(lamb)
                Ind=length(lamb)-1;
                
            end
            lref=1i*(imag(lamb(1,Ind-1)):imag(lamb(1,Ind+1)-lamb(1,Ind-1))/200:imag(lamb(1,Ind+1)));
            
SR=SpectralRadius(W,lref,'real');
r0=min(SR);

res.lam=lambda;
res.R=r;
res.R0=r0;

end