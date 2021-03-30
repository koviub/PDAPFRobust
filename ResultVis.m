function res=ResultVis(fn,a,q,rend)

load(fn);

tt=res.tau(1,:);
tk=res.tau;
lk=res.Length;
sstab=(1-a)/a;
if isempty(res.Radii)
    pl=log(res.Area);
else
    pl=res.Radii;
end


step=50;
[ttk,lt]=meshgrid(0.05:.95/step:1,0:10/step:10);
f=interp2(tk,lk,pl,ttk,lt,'cubic');
figure()

hold on;
surface(ttk(2:end-2,2:end-2),f(2:end-2,2:end-2),lt(2:end-2,2:end-2))

% [~,h]=contour(ttk(2:end-2,2:end-2),f(2:end-2,2:end-2),lt(2:end-2,2:end-2),[.1 .3 .5 1 2 3]);
% contour(ttk(2:end-2,2:end-2),f(2:end-2,2:end-2),lt(2:end-2,2:end-2),[0.3 .5],'r')
zlim([0 14])
caxis([0 3])
ylim([0 1])
xlim([.05 .5])
h.ShowText='on';
set(gca,'fontsize',18)
ylabel('$r^{\bf\Delta}$','interpreter','latex')
xlabel('$\tau$[s]','interpreter','latex')
grid on; box on;


switch nargin
    case 1
        a=0;
        q=1;
        rend=max(max(f));
    case 2
        q=1;
        rend=max(max(f));
    case 3
        rend=max(max(f));
end

if q~=0
    ft=fittype('a1*x.^2');
else
    ft=fittype('a1*x');
end

exc=1;
a0=(3/(4*(1+a))*9.81);
h=1;
valEps(1)=0;
valA1(1)=a0;
valR(1)=0;

figure()
hold on;grid on;
cmap=colormap;
%(length(cmap)-1)
for radi=0.001:(rend-.001)/50:rend
    h=h+1;
    fi=f;
    fi(fi>=radi)=0;
    [~,ind]=max(fi);
    for kk=1:size(ttk,1)
        valL(kk,h)=lt(ind(kk),kk);
        valT(kk,h)=ttk(ind(kk),kk);
    end
    Tau=zeros(length(valT(:,h)),1);
        Tau(1:end)=valT(:,h);
        LL=zeros(length(valL(:,h)),1);
        LL(1:end)=valL(:,h);
    exc=(LL(:)>=6)|(LL(:)<=1.1);%|(LL(:)<=.36)|(Tau(:)>.425)|(LL(:)>5*a0*Tau(:).^2);%|(LL(:)<1/4*a0*Tau(:).^2);%
%      exc=(LL(:)>=8)|(LL(:)<=.06)|(LL(:)>5*a0*Tau(:).^2)|(LL(:)<1/4*a0*Tau(:).^2);
    % if wrong fit check exclude data
    [fitted,gof,~]=fit(Tau,LL,ft,'Exclude',exc,'startpoint',[1]);
    valA1(h)=fitted.a1;
    valEps(h)=radi;
%         hold on
%         plot(fitted,Tau,LL,exc)
        LLC=fitted.a1*Tau(:).^2;
        plot(Tau,LLC,'color',cmap(h-1,:))
        legend('off')
        valR(h)=gof.rsquare;
end

index=find(sstab>valEps);
valEps(index(end):end)=sstab;
% contour(tk,lk,pl,0.1:(rend-.1)/10:rend);
% plot(f,[t l],rad,'Exclude',indi)
% [cc,hh]=contour(tk,lk,pl,[1 1]*.5);
% plot(tt,3/(4*(1+a))*9.81*tt.^2,'r')
plot(tt,3./(4*(1+a)).*9.81.*tt.^2,'color',cmap(1,:))
ylim([0 10])
xlabel('$\tau$','interpreter','latex')
ylabel('$L_{\rm crit}$','interpreter','latex')

[tt,ep]=meshgrid(Tau(1:length(Tau)/49:end),valEps);
ep(ep==Inf)=1;
LL=[];
for i=1:size(tt,2)
    for j=1:size(ep,1)
        LL(j,i)=valA1(j)*tt(1,i)^2;
    end
end
figure();
hold on
contour(tt(2:end-2,2:end-2),ep(2:end-2,2:end-2),LL(2:end-2,2:end-2),[0.12 1 2 3 4 5 6 7])
contour(tt(2:end-2,2:end-2),ep(2:end-2,2:end-2),LL(2:end-2,2:end-2),[.3 .5],'r')
zlim([0 14])
caxis([0 7])
ylim([0 2])
xlim([0.1 1])
set(gca,'fontsize',18)
ylabel('$r^{\bf\Delta}$','interpreter','latex')
xlabel('$\tau$[s]','interpreter','latex')
% 
% plot([.2 .2],[0 .6],'k')
% plot([.0 .2],[0.041 .041],'k')
% plot([.25 .25],[0 .6],'k')
% plot([.0 .25],[0.0024168 0.0024168],'k')
plot([.2 .2],[0 .6],'k')
plot([.25 .25],[0 .6],'k')
box on
grid on

figure();
contour(tt,ep,LL)


figure('position',[10 10 400 400],'Name',[fn(1:end-4) '_a1'])
hold on;grid on;

plot(valEps,valA1./a0,'r')
ylabel('$a_1/a_0$','interpreter','latex')
xlabel('$\epsilon$','interpreter','latex')
res.f=f;res.tt=ttk;res.lt=lt;res.a1=valA1./a0;res.eps=valEps;

% if strcmp(fn(1:6),'system')
%     ft1=fittype('a*x+1');
% else
%     ft1=fittype('exp(a*x^2)');
% end
% rf=fit(res.eps',res.a1',ft1);%,'Exclude',res.eps>.8
% plot(rf,'b')
xlim([0 max(valA1./a0)])
% res.rf=rf;
end