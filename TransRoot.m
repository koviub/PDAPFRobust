function [ r ] = TransRoot( Dcar,lambda, order, h, R, draw )
%TransRoot Roots of transcendental equations.
%   TransRoot(D, s, N, X1, R, B) calculates N roots of a transcendental 
%   equation D(s)=0, within a closed region in the complex plane, 
%   with center X1 and radius R.
%
%   See reference: http://dx.doi.org/10.1155/2015/523043
warning('off')
if nargin<6
    draw=false;
end

samp=2^8;
T=1/samp;
L=2*pi*T;
th=(0:(samp-1))*L;

fi=0:pi/24:2*pi;
x=real(h)+R*cos(fi);
y=imag(h)+R*sin(fi);

ll=R*exp(1i*th)+h;

F=subs(Dcar,lambda,ll);
G=double(1./F);

gg=fourier(G);
g=gg(2:2*order+1);

for ii=1:order
    H(1:order,ii)=g(ii:(ii-1+order));
end

B(:,1)=-g(order+1:end);

c=linsolve(H,B);
c(end+1)=1;

c=flip(c);
r=R*roots(c)+h;

if draw
    hold on
    plot(r,'r*')
    plot(x,y,'k:')
    
% contours of real and imaginary parts are zero
    xx=linspace(real(h)-R,real(h)+R,30);
    yy=linspace(imag(h)-R,imag(h)+R,30);
    [ga,o]=meshgrid(xx,yy);
    DD=double(subs(Dcar,lambda,ga+1i*o));
    RE=real(DD);
    IM=imag(DD);
    
    contour(ga,o,RE,[0 0],'r:')
    contour(ga,o,IM,[0 0],'b:')
end
end