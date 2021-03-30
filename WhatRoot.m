function varargout = WhatRoot(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WhatRoot_OpeningFcn, ...
    'gui_OutputFcn',  @WhatRoot_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before WhatRoot is made visible.
function WhatRoot_OpeningFcn(hObject, ~, handles, varargin)
% Choose default command line output for WhatRoot
global input zoom anal TZM pseudo sol

input=handles;
zoom=get(handles.radiobutton1,'Value');
TZM = get(handles.radiobutton3,'Value');
anal = get(handles.radiobutton2,'Value');
pseudo = get(handles.radiobutton4,'Value');

handles.output = hObject;
set(handles.axes1,'box','on','xgrid','on','ygrid','on')
set(handles.axes2,'box','on','xgrid','on','ygrid','on')
set(handles.pushbutton3,'String','Delay margin');
set(handles.text11,'backgroundcolor','white')
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WhatRoot wait for user response (see UIRESUME)
% uiwait(handles.figure1);

sol = AnalyticSolu(handles);
end
% --- Outputs from this function are returned to the command line.
function varargout = WhatRoot_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input zoom anal sol TZM pseudo

input=handles;
zoom=get(handles.radiobutton1,'Value');
TZM = get(handles.radiobutton3,'Value');
anal = get(handles.radiobutton2,'Value');
pseudo = get(handles.radiobutton4,'Value');
var1=str2num(get(handles.edit2,'String'));
var2=str2num(get(handles.edit3,'String'));
var3=str2double(get(handles.edit4,'String'));
var4=str2double(get(handles.edit5,'String'));

addpath('C:\Users\Kovacs Balazs\Desktop\MDBM-Matlab-master\code_folder')


var=str2num(get(handles.edit6,'String'));

par.params=str2num(get(handles.edit6,'String'));
param=par.params;
axes(handles.axes1)
cla
box on;hold on;

if get(handles.radiobutton6,'value')%
    
    if var1(2)>sol.Var01&&var2(2)>sol.Var02&&var1(1)<sol.Var01&&var2(1)<sol.Var02
        xlim([var1(1),var1(2)])
        ylim([var2(1),var2(2)])
    else
        x1=sol.Var01*abs(sign(sol.Var01)-.05);y1=sol.Var02*abs(sign(sol.Var02)-.05);
        x2=sol.Var01*5;y2=abs(sol.Var02)*05;
        str=sprintf('[%f,%f]',x1,x2);
        set(handles.edit2,'String',str);
        str=sprintf('[%f,%f]',y1,y2);
        set(handles.edit3,'String',str);
        
        xlim([x1,x2])
        ylim([y1,y2])
    end
else
    if .5*3*9.81/4<param(1)/param(2)^2
        x1=sol.Var01*.95;
        eq0=@(x)double(-sol.VarC2(x));
        om=0.1:.1:110;
        [dk0,~]=findpeaks(eq0(om));
        dk0=-dk0;
        eq1=@(x)double(sol.VarC1(x));
        [pk1,ind]=findpeaks(eq1(om));
        if sol.Var02<dk0(1)
            y1=sol.Var02*abs(sign(sol.Var02)-.05);
        else
            y1=dk0(1)*abs(sign(dk0(1))-.05);
        end
        eq2=@(x)double([sol.VarC1(x(1))-sol.VarC1(x(2));sol.VarC2(x(1))-sol.VarC2(x(2))]);
        om2=fsolve(eq2,om(ind(1))*[1;3]);
        pk2=double(sol.VarC1(om2(1)));
        if pk1(1)<4*pk2
            x2=pk1(1)*1.01;
        else
            x2=pk2*1.01;
        end
        eq3=@(x)double(sol.VarC1(x)-sol.Var01);
        om3=fsolve(eq3,1.5*om2(2));
        y2=abs(double(sol.VarC2(om3))*1.1);
        str=sprintf('[%f,%f]',x1,x2);
        set(handles.edit2,'String',str);
        str=sprintf('[%f,%f]',y1,y2);
        set(handles.edit3,'String',str);
        
        xlim([x1,x2])
        ylim([y1,y2])
    else
        om=0.1:.1:110;
        x1=sol.Var01*abs(sign(sol.Var01)-.2);y1=sol.Var02*abs(sign(sol.Var02)-.2);
        eq1=@(x)double(sol.VarRef1(x));
        [pk,ind]=findpeaks(eq1(om));
        eq2=@(x)double(sol.VarRef1(x)-sol.VarRef01);
        dk=fzero(eq2,om(ind(1)));
        x2=pk(1)*1.01;y2=double(sol.VarRef2(dk))*1.01;
        str=sprintf('[%f,%f]',x1,x2);
        set(handles.edit2,'String',str);
        str=sprintf('[%f,%f]',y1,y2);
        set(handles.edit3,'String',str);
        
        xlim(sort([x1 x2]))
        ylim(sort([y1 y2]))
    end
end


var1=str2num(get(handles.edit2,'String'));
var2=str2num(get(handles.edit3,'String'));

ax=[];
ax(1).val=linspace(var1(1),var1(2),var4);  % var1
ax(2).val=linspace(var2(1),var2(2),var4);  % var2
ax(3).val=linspace(0,var3,3*var4);  % om
       
if ~anal
    v1=sol.VarR1(ax(2).val).*ones(1,length(ax(2).val));
    v2=sol.VarR2(ax(2).val);
    v3=double([sol.Var01 sol.VarC1(ax(3).val(2:end))]);
    v4=double([sol.Var02 sol.VarC2(ax(3).val(2:end))]);
    
    
    tol = 100; % distance > tol indicates discontinuity
    dl = diff([v3;v4],1,2); % look up what that command does if you don't know it
    euler_dist = sqrt((dl(1,:)+dl(2,:)).^2); % distance between data points
    jumpind = [0 euler_dist>tol]; % now if jumpind(i) = true, we know that the
    %   point [lat(i) lon(i)] is the first after a jump
    blocks = cumsum(jumpind); % points that belong to the same continuous part
    % have the same value in blocks
    % Now just loop over the continuous blocks to draw a separate line for each one
    for i=0:blocks(end)
        plot(v3(blocks==i),v4(blocks==i),'b');
        hold on;
    end
    plot(v1,v2,'r')
else
    par.interpolationorder=3;
    Niteration=4;
    bound_function_name='Model';
    mdbm_sol=mdbm(ax,bound_function_name,Niteration,[],par);
    plot_mdbm(mdbm_sol);
    view([0,0,1])
end
xlabel('var1')
ylabel('var2')

end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(~, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pseudo
pseudo=get(handles.radiobutton4,'Value');
set(handles.text11,'String','')

param=str2num(get(handles.edit6,'String'));
var=str2num(get(handles.edit7,'String'));

axes(handles.axes1)
box on;hold on;
plot(var(1),var(2),'ro','markerfacecolor','r')

axes(handles.axes2)
cla
box on;hold on;

x0=str2double(get(handles.edit8,'String'));
Rad=str2double(get(handles.edit10,'String'));
n_r=str2double(get(handles.edit11,'String'));

F=str2func(get(handles.edit1,'String'));
syms s
Dcar=F(s,var,param);
r=TransRoot(Dcar,s,n_r,x0,Rad,true);

plot([0 0],[imag(x0)-Rad imag(x0)+Rad],'k--','linewidth',1.5)
for i=1:length(r)
    if real(r(i))>0
        cb='r';
    else
        cb='g';
    end
    if abs(double(subs(Dcar,r(i))))<1e-2
        str=[get(handles.text11,'String') sprintf('%.2f+i*%.2f, ',real(r(i)),imag(r(i)))];
        set(handles.text11,'String',str)
        plot(r(i),'ro','markersize',6,'markerfacecolor',cb)
        drawnow();
    end
    
end

if pseudo
    res=Robustness(handles,true);
    Re=real(res.lam);
    Im=imag(res.lam);
    R=res.R;
    
    contour(Re,Im,R,0:.2:1)
    contour(Re,Im,R,[1 1].*res.R0,'k')
    
    str=[get(handles.text11,'String') sprintf('r=%.2f, ',res.R0)];
    set(handles.text11,'String',str)
end
xlim([real(x0)-1.1*Rad real(x0)+1.1*Rad])
ylim([imag(x0)-1.1*Rad imag(x0)+1.1*Rad])
xlabel('Re(s)')
ylabel('Im(s)')

end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pseudo rmax
if ~pseudo
    res = MarginofStab(handles);
    var=str2num(get(handles.edit6,'String'));
    rmax=0;
    var(2)=res.x;
    str='[';
    for i=1:length(var)
        str=[str sprintf('%f,',var(i))];
    end
    str(end)=']';
    set(handles.edit6,'String',str);
    edit6_Callback(hObject, eventdata, handles);
    pushbutton1_Callback(hObject, eventdata, handles);
    
else
    res = SearchMax(handles);
    rmax=res.Rmax;
end
x=res.var1;y=res.var2;
str=sprintf('[%f,%f]',x,y);
set(handles.edit7,'String',str);

pushbutton2_Callback(hObject, eventdata, handles);
% plot(x,y,'ro','markerfacecolor','r')
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sol rmax pseudo
par=str2num(get(handles.edit6,'String'));
par1=str2num(get(handles.edit12,'String'));
par2=str2num(get(handles.edit13,'String'));
[x,y]=meshgrid(par1,par2);
rM=zeros(size(x));
nmax=length(par2)*length(par1);
cmap=parula(256);

if get(handles.radiobutton5,'Value')
    pushbutton1_Callback(hObject, eventdata, handles);
    pushbutton3_Callback(hObject, eventdata, handles);
else
    %     fig=figure();
    %     hold on;grid on;box on;
    %     pl1=plot(x(:),y(:),'*','color',cmap(1,:));
    str=sprintf('[%f,%f]',sol.Var01,sol.Var02);
    set(handles.edit7,'String',str)
    if pseudo
        i=2;n=1;
        %     wh=waitbar(i,sprintf('Wait...,%d / %d',n,nmax));
        str=sprintf('%.2f, ',n/nmax);
        set(handles.text16,'String',str)
        for p1=par1
            j=9;
            for p2=par2
                % cycle through par1 and par2
                set(handles.edit6,'String',sprintf('[%f,%f,%f,%f]',p1,p2,par(3),par(4)));
                sol=AnalyticSolu(handles);
                pushbutton1_Callback(hObject, eventdata, handles);
                par1=str2num(get(handles.edit12,'String'));
                par2=str2num(get(handles.edit13,'String'));
                str=sprintf('[%f,%f]',sol.Var01,sol.Var02);
                set(handles.edit7,'String',str);
                pushbutton3_Callback(hObject, eventdata, handles);
                rM(j,i)=rmax;
                %             c_index=floor(rmax*255)+1;
                %             set(pl1,'xdata',x(i,j),'ydata',y(i,j),'color',cmap(c_index,:))
                %             waitbar(n/nmax,wh,sprintf('Wait...,%d / %d',n,nmax))
                
                % save data
                filename=lower(sprintf('PF_%.2f_%.2f',par(3),par(4)));
                saveas(handles.figure1,[sprintf('Fig_%d_%d',i,j) filename '.png'])
                str=sprintf('%.2f, ',n/nmax);
                set(handles.text16,'String',str)
                j=j+1;n=n+1;
            end
            i=i+1;
        end
        res.x=x;res.y=y;res.z=rM;
        
    else
        % calculate margin on plane par1, par2 for fixed par(3), par(4)
    end
    %     close(wh)
    filename=lower(sprintf('PF_%.2f_%.2f',par(3),par(4)));
    save(['Data_' filename '.mat'],'res')
end
end

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zoom

axes(handles.axes1)
box on;hold on;
if zoom
    a=getrect(handles.axes1);
    if a(3)>0&&a(4)>0
        x1=a(1);y1=a(2);
        x2=x1+a(3);y2=y1+a(4);
        str=sprintf('[%f,%f]',x1,x2);
        set(handles.edit2,'String',str);
        str=sprintf('[%f,%f]',y1,y2);
        set(handles.edit3,'String',str);
        
        xlim([x1,x2])
        ylim([y1,y2])
    else
        
        x1=-5;y1=-5;
        x2=25;y2=25;
        str=sprintf('[%f,%f]',x1,x2);
        set(handles.edit2,'String',str);
        str=sprintf('[%f,%f]',y1,y2);
        set(handles.edit3,'String',str);
        xlim([-5,25])
        ylim([-5,25])
    end
    
else
    [x,y]=getpts(handles.axes1);
    str=sprintf('[%f,%f]',x,y);
    set(handles.edit7,'String',str);
    
    plot(x,y,'ro','markerfacecolor','r')
    pushbutton2_Callback(hObject, eventdata, handles);
end
end


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2)
box on;hold on;
a=getrect(handles.axes2);
if a(3)>0&&a(4)>0
    x1=a(1);y1=a(2);
    x2=x1+a(3);y2=y1+a(4);
    xlim([x1,x2])
    ylim([y1,y2])
else
    axis auto
end
end

function edit1_Callback(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% str1='@(x,var,par)x^2+2*par(2)*par(1)*x+par(1)^2+var(2)-var(2)*exp(-x*60/var(1))';
% str1='@(x,var,par)x^2-3*9.81/2/par(1)+var(1)*exp(-par(2)*x)+var(2)*x*exp(-par(2)*x)+par(4)*x^2*exp(-par(2)*x)';
str1='@(x,var,par)par(4)*x^2*exp(-par(2)*x)-3*9.81/2/par(1) + x*(var(2) + x) + var(1) + (2*(-(sqrt(3*9.81/2/par(1))*(var(2)*x + var(1))*cosh(sqrt(3*9.81/2/par(1))*(1 + par(3))*par(2))) - (3*9.81/2/par(1)*var(2) + x*var(1))*sinh(sqrt(3*9.81/2/par(1))*(1 + par(3))*par(2)))*sinh((x*(par(2) - (1 + par(3))*par(2)))/2))/(sqrt(3*9.81/2/par(1))*exp((x*(par(2) + (1 + par(3))*par(2)))/2))';
set(hObject,'String',str1);

end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
end
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str1='[-5,25]';
set(hObject,'String',str1);
end

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

end
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str1='[-5,25]';
set(hObject,'String',str1);
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
end
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str1='50';
set(hObject,'String',str1);
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
end
% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
str1='31';
set(hObject,'String',str1);
end

function edit6_Callback(hObject, ~, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
global sol

sol = AnalyticSolu(handles);
end
% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% str1='[3*9.81/2/2,0.1]';
str1='[1,0.2,-1,.0]';
set(hObject,'String',str1);
end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

end
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

pushbutton2_Callback(hObject, eventdata, handles);
end
% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

pushbutton2_Callback(hObject, eventdata, handles);
end
% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


end

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

pushbutton2_Callback(hObject, eventdata, handles);
end
% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

pushbutton2_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','0:.8:8');
end

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','0:.1:1');
end

function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double

end

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','[1,1,0]')
end

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
global zoom

zoom = get(hObject,'Value');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
global anal

anal = get(hObject,'Value');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
global TZM

TZM = get(hObject,'Value');
end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
global pseudo

pseudo = get(hObject,'Value');

if pseudo
    set(handles.pushbutton3,'String','Robust max');
else
    set(handles.pushbutton3,'String','Delay margin');
end

end

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
end


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
end
