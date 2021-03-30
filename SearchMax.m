function res = SearchMax(input)
syms var1 var2 real
syms s

var=[var1 var2];
nm=str2double(get(input.edit5,'String'));
cmap=parula(256);

param=str2num(get(input.edit6,'String'));
weights=str2num(get(input.edit16,'String'));
vari=str2num(get(input.edit7,'String'));
x0=str2double(get(input.edit8,'String'));                          % origin on the complex plane
Rad=str2double(get(input.edit10,'String'));
n_r=str2double(get(input.edit11,'String'));
F=str2func(get(input.edit1,'String'));

vari1=str2num(get(input.edit2,'String'));
vari2=str2num(get(input.edit3,'String'));
Vv1=linspace(vari1(1),vari1(2),nm);
Vv2=linspace(vari2(1),vari2(2),nm);
[V1,V2]=meshgrid(Vv1,Vv2);
r=zeros(size(V1));
ix0=(3*vari1(1)+vari1(2))/4;%sum(vari1)/2;% vari(1)
iy0=(3*vari2(1)+vari2(2))/4;%sum(vari2)/2;% vari(2)
i0=floor((ix0-vari1(1))/(vari1(2)-vari1(1))*nm)+1;j0=floor((iy0-vari2(1))/(vari2(2)-vari2(1))*nm)+1; % first point from picked / transform values to indeces!
i1=[0 1 1 1 0 -1 -1 -1 0 1 2 2 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1 0 3 0 -3];
j1=[1 1 0 -1 -1 -1 0 1 2 2 2 1 0 -1 -2 -2 -2 -2 -2 -1 0 1 2 2 3 0 -3 0];
% i1=[0 1 1 1 0 -1 -1 -1 0 2 0 -2];
% j1=[1 1 0 -1 -1 -1 0 1 2 0 -2 0];
if i0<1
    i0=1;
elseif i0>nm
    i0=nm;
end

if j0<1
    j0=1;
elseif j0>nm
    j0=nm;
end
inext=i0;jnext=j0;
% caulate first
% structured perturbation on var1 var2
% covar=coeffs(F(s,var,param),var);
Dc=@(x,v)F(x,v,param);
% w1=@(x)double(subs([weights(1)*vari(1)*covar(3);weights(2)*vari(2)*covar(2);weights(3)*param(4)*s^2*exp(-param(2)*s)],'s',x));
% w1=@(x)double([vari(1);vari(2)*x;param(4)*x^2].*exp(-param(2)*x));
w1=@(x,v)double([weights(1)*v(1)*(1 - (exp(-(x*(param(2) + param(2)*(param(3) + 1)))/2)*sinh((x*(param(2) - param(2)*(param(3) + 1)))/2)*(2*x*sinh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1)) + 2*(3/2*9.81/param(1))^(1/2)*cosh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1))))/(3/2*9.81/param(1))^(1/2));...
    weights(2)*v(2)*(x - (exp(-(x*(param(2) + param(2)*(param(3) + 1)))/2)*sinh((x*(param(2) - param(2)*(param(3) + 1)))/2)*(2*(3/2*9.81/param(1))*sinh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1)) + 2*(3/2*9.81/param(1))^(1/2)*x*cosh((3/2*9.81/param(1))^(1/2)*param(2)*(param(3) + 1))))/(3/2*9.81/param(1))^(1/2));...
    weights(3)*param(4)*x^2.*exp(-param(2)*x)]);% pf weight function

W=@(x)double(w1(x,vari)./Dc(x,vari));

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

axis(input.axes1);
c_index=floor(r0/20*255)+1;
plot(vari(1),vari(2),'color',cmap(c_index,:),'Marker','o','MarkerFaceColor',cmap(c_index,:),'markersize',2)
drawnow();

% step in the directon of the largest difference while maximum is found
% find maximum
while true
    for f=1:length(i1)
        indi=i0+i1(f);
        indj=j0+j1(f);
        
        % bound to viewport
        if indi<1
            indi=1;
        elseif indi>nm
            indi=nm;
        end
        
        if indj<1
            indj=1;
        elseif indj>nm
            indj=nm;
        end
        
        if r(indj,indi)==0
            fprintf('*');
            vari=[V1(indj,indi),V2(indj,indi)];
            
            W=@(x)w1(x,vari)./Dc(x,vari);
            
            root=TransRoot(Dc(s,vari),s,n_r,x0,Rad,false);
            for i=1:length(root)
                if real(root(i))>0
                    cb(i)=1;
                else
                    cb(i)=0;
                end
                if abs(double(Dc(root(i),vari)))>1e-2
                    cb(i)=0;
                end
            end
            
            if sum(cb(:))==0
                [~,Ind]=min(SpectralRadius(W,lamb,'real'));
                if Ind==1
                    Ind=2;
                elseif Ind==length(lamb)
                    Ind=length(lamb)-1;
                    
                end
                lref=1i*(imag(lamb(1,Ind-1)):imag(lamb(1,Ind+1)-lamb(1,Ind-1))/1000:imag(lamb(1,Ind+1)));
                r(indj,indi)=min(SpectralRadius(W,lref,'real'));
            else
                r(indj,indi)=0;
            end
        end
        
        %plot here on axis1
        axis(input.axes1);
        c_index=floor(r(indj,indi)*255)+1; %/max(r(:))
        plot(vari(1),vari(2),'color',cmap(c_index,:),'Marker','o','MarkerFaceColor',cmap(c_index,:),'markersize',2)
        drawnow();
        
        if r(j0,i0)<r(indj,indi)
            inext=indi;jnext=indj;
        end
    end
    fprintf('\n');
    if inext==i0&&jnext==j0
        rmax=r(j0,i0);
        var_max=[V1(j0,i0),V2(j0,i0)];
        fprintf('found max');
        break;
    end
    
    i0=inext;j0=jnext;
end

res.Rstart=r0;
res.Rpath=r;
res.Rmax=rmax;
res.var1=var_max(1);
res.var2=var_max(2);

end