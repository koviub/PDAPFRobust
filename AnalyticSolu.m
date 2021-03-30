function sol = AnalyticSolu(input)
syms var1 var2 om et real

var=[var1 var2];
lam=1i*om;

pp=str2num(get(input.edit6,'String'));
param=[pp(1) pp(2) et pp(4)];
F=str2func(get(input.edit1,'String'));
Dcar0=F(0,var,param);
Dcar1=F(lam,var,param);

Eq(1)=Dcar0==0;
Eq(2)=real(Dcar1)==0;
Eq(3)=imag(Dcar1)==0;

sol0=solve(Eq(1),var1);

VarR1=@(var2)double(sol0);
VarR2=@(var2)double(var2);

sol1=solve(Eq(2:3),var);

VarC1=@(s)subs(sol1.var1,{'om','et'},{s,pp(3)});
VarC2=@(s)subs(sol1.var2,{'om','et'},{s,pp(3)});
VarRef1=@(s)subs(sol1.var1,{'om','et'},{s,-1});
VarRef2=@(s)subs(sol1.var2,{'om','et'},{s,-1});


VarRef01=double(limit(VarRef1(om),om,0));
VarRef02=double(limit(VarRef2(om),om,0));
Var01=double(limit(VarC1(om),om,0));
Var02=double(limit(VarC2(om),om,0));

sol=struct('VarR1',VarR1,'VarR2',VarR2,...
    'VarRef1',VarRef1,'VarRef2',VarRef2,'VarRef01',VarRef01,'VarRef02',VarRef02,...
    'VarC1',VarC1,'VarC2',VarC2,'Var01',Var01,'Var02',Var02);

end
