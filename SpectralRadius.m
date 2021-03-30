function rR = SpectralRadius(W,lambda,str)
%% Value of spectral radius on each point of a complex region
% input: W      -- Resolvent matrix
%        lambda -- point on complex plane 
%        str    -- perturbation type
%

if nargin<3
    str='complex';
end

ss=size(lambda);
rR=zeros(ss(1),ss(2));
switch lower(str)
    case 'complex'
        for i=1:size(lambda,1)
            for j=1:size(lambda,2)
                v=norm(W(lambda(i,j)));
                rR(i,j)=1/v;
            end
        end
        
    case 'real'
        for i=1:size(lambda,1)
            for j=1:size(lambda,2)
                v=muR(W(lambda(i,j)));
                rR(i,j)=1/v;
            end
        end
end

end

function v=muR(W)
%%
val=zeros(10,1);
i=1;

for gamma=.0001:.1:1
Mx=[real(W),-gamma*imag(W);...
    1/gamma*imag(W),real(W)];

%% second largest singular value
test1=isnan(Mx);
test2=Mx==Inf;
if sum(test1(:))==0&&sum(test2(:))==0
    vec=SingularValues(Mx);
    % vec= svds(Mx,2); % slow calculates singular value directions also...
else
    vec=zeros(2,1);
end
val(i)=vec(2);
i=i+1;
end

v=min(val);
end

function sing_values=SingularValues(Mx)

if (ismatrix(Mx)&&diff(size(Mx))) %matrix is not square
% sing_values=svds(Mx,2,'largest');

A=Mx.'*Mx;
sing_values=sqrt(sort(abs(eig(A)),'descend'));

else % for square matrix 
sing_values=sort(abs(eig(Mx)),'descend');
end

end