function res = MarginofStab(input)
syms var1 var2 om tau real
syms s
global sol TZM

var=[var1 var2];
lam=1i*om;

p=str2num(get(input.edit6,'String'));
param=[p(1) tau p(3) p(4)];
F=str2func(get(input.edit1,'String'));
F0=eval(['@(s,var1,var2,tau)' char(F(s,var,param))]);
F1=eval(['@(s,var1,var2,tau)' char(diff(F0(s,var1,var2,tau),s))]);
F2=eval(['@(s,var1,var2,tau)' char(diff(F1(s,var1,var2,tau),s))]);

Eq(1,1)=F0(0,var1,var2,tau);
Eq(2,1)=F1(0,var1,var2,tau);
Eq(3,1)=F2(0,var1,var2,tau);
Eq(4,1)=real(F0(lam,var1,var2,tau))==0;
Eq(5,1)=imag(F0(lam,var1,var2,tau))==0;

sol1=solve(Eq(4:5),var);

VarC1=@(s,t)subs(sol1.var1,[om,tau],[s,t]);
VarC2=@(s,t)subs(sol1.var2,[om,tau],[s,t]);

Var01=@(tau)limit(VarC1(om,tau),om,0);
Var02=@(tau)limit(VarC2(om,tau),om,0);

if TZM
% derivative at om~0
% F1=eval(['@(s,tau)' char(diff(VarC1(s,tau),s))]);
% F2=eval(['@(s,tau)' char(diff(VarC2(s,tau),s))]);
% 
% slp=@(tau)limit(F1(om,tau)/F2(om,tau),om,0);
% eq=@(x)double(slp(x));
% x0=p(2);
% x1=fsolve(eq,x0);
% tres=x1;

eq=@(x)double([F0(0,x(1),x(2),x(3));...
    F1(0,x(1),x(2),x(3));...
    F2(0,x(1),x(2),x(3))]);
x0=double([Var01(p(2)),Var02(p(2)),p(2)]);

x1=fsolve(eq,x0);

res.var1=x1(1);
res.var2=x1(2);
tres=x1(3);
end

% ComplexDZM
eq=@(x)double([Var01(x(2))-VarC1(x(1),x(2));...
            Var02(x(2))-VarC2(x(1),x(2))]);
x0=[10;p(2)];
if ~TZM
x2=fsolve(eq,x0);
tres=x2(2);
end
res.x=tres;

res.var1=Var01(res.x);
res.var2=Var02(res.x);

end
