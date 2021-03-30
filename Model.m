function Eq = Model(ax,par)
global input
NumEq=2;
NumVar=size(ax,2);
Eq=zeros(NumEq,NumVar);
for kax=1:NumVar
    p1=ax(1,kax);
    p2=ax(2,kax);
    om=ax(3,kax);
    lam=1i*om;
    var=[p1,p2];
    
    param=par.params;
    
    F=str2func(get(input.edit1,'String'));
    Dcar=F(lam,var,param);
    
    Eq(1,kax)=real(Dcar);
    Eq(2,kax)=imag(Dcar);
end
end
