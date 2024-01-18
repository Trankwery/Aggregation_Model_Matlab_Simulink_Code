function varargout = RenderDzons_3D(varargin)
% The function  RenderDzons_3D() allows the visualization 
% of the result of the simulation. 
% 
% The input parameters of the function are as follows: 
%  RenderDzons_3D(Pout,Rd,Rinkl)
%    Pout  - the result of the simulations
%    Rd    - the radius of the droplet
%    Rinkl - the radius of inclusion
% 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RenderDzons_3D_OpeningFcn, ...
                   'gui_OutputFcn',  @RenderDzons_3D_OutputFcn, ...
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


% --- Executes just before RenderDzons_3D is made visible.
function RenderDzons_3D_OpeningFcn(hObject, eventdata, handles, varargin)
% varargin   command line arguments to RenderDzons_3D (see VARARGIN)

% Choose default command line output for RenderDzons_3D
handles.output = hObject;
handles.poz = varargin{1}.signals.values;
handles.time = varargin{1}.time;
handles.Rk = varargin{2};
handles.Rinkl = varargin{3};
handles.length = size(handles.poz,3);
set(handles.slider1,'max',handles.length,'value',1,'min',1,'sliderstep',[1/handles.length,10/handles.length]);
% Update handles structure
cameratoolbar;
if ~isscalar( handles.Rk )
    [C I] =  min( abs(handles.Rk(:,1) - handles.time( 1 )) );
    Rk = handles.Rk(I,2);
%     set(handles.teTime,'String',['Carent time : ' num2str(handles.time( 1 ) ) ] );
else
    Rk = handles.Rk;
end
plot3([-1.2*Rk 1.2*Rk],[0 0],[0 0]);
hold on;
grid on;
plot3([0 0],[-1.2*Rk 1.2*Rk],[0 0]);
plot3([0 0],[0 0],[-1.2*Rk 1.2*Rk]);
guidata(hObject, handles);
% UIWAIT makes RenderDzons_3D wait for user response (see UIRESUME)
function varargout = RenderDzons_3D_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
nom = floor( get(handles.slider1,'value') );
R = handles.Rinkl; 

if ~isscalar( handles.Rk )
%     [C I] =  min( abs(handles.Rk(:,1) - handles.time( nom )) );
    rd = handles.Rk(nom,2);
%     set(handles.teTime,'String',['Carent time : ' num2str(handles.time( nom ) ) ] );
else
    rd = handles.Rk;
end

        xx = rd * sin(0:0.01:2*pi);
        yy = rd * cos(0:0.01:2*pi);
        zz = zeros(1,length(xx));
        try
            set(handles.kr(1),'xdata',xx,'ydata',yy,'zdata',zz);
            set(handles.kr(2),'xdata',xx,'ydata',zz,'zdata',yy);
            set(handles.kr(3),'xdata',zz,'ydata',xx,'zdata',yy);
        catch
            handles.kr(1) = plot3(xx,yy,zz);
            handles.kr(2) = plot3(xx,zz,yy);
            handles.kr(3) = plot3(zz,xx,yy);
            guidata(hObject, handles);
        end
            u = (0:0.1*pi:2*pi)';
            v = [0:0.1*pi:2*pi];
        for ii = 1 : size(handles.poz,1)
            X1 = R*sin(u)*cos(v) + handles.poz(ii,1,nom);
            Y1 = R*sin(u)*sin(v) + handles.poz(ii,2,nom);
            Z1 = R*cos(u)*ones(size(v)) + handles.poz(ii,3,nom);
            try
            set(handles.hS(ii),'xdata',X1,'ydata',Y1,'zdata',Z1);
            catch
                handles.hS(ii) = mesh(X1,Y1,Z1);
                guidata(hObject, handles);
            end
        end
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function meRotate_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Zoom_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function meTools_Callback(hObject, eventdata, handles)

