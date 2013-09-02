function varargout = rrgmres_demo(varargin)
% RRGMRES_DEMO M-file for rrgmres_demo.fig
%       RRGMRES_DEMO, by itself, creates a new RRGMRES_DEMO or raises the 
%       existing instance.
%
%       H = RRGMRES_DEMO returns the handle to a new RRGMRES_DEMO or the 
%       handle to the existing instance.
%
%       RRGMRES_DEMO('CALLBACK',hObject,eventData,handles,...) calls the 
%       local function named CALLBACK in RRGMRES_DEMO.M with the given 
%       input arguments.
%
%       RRGMRES_DEMO('Property','Value',...) creates a new RRGMRES_DEMO or 
%       raises the existing instance.  Starting from the left, property 
%       value pairs are applied to the GUI before rrgmres_demo_OpeningFcn 
%       gets called.  An unrecognized property name or invalid value makes 
%       property application stop.  All inputs are passed to 
%       rrgmres_demo_OpeningFcn via varargin.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RRGMRES_DEMO

% Last Modified by GUIDE v2.5 14-May-2011 11:59:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sym_rrgmres_demo_OpeningFcn, ...
                   'gui_OutputFcn',  @sym_rrgmres_demo_OutputFcn, ...
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
% --- Executes just before RRGMRES_DEMO is made visible.
function sym_rrgmres_demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RRGMRES_DEMO (see VARARGIN)

% Choose default command line output for RRGMRES_DEMO
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using RRGMRES_DEMO.
if strcmp(get(hObject,'Visible'),'off')
    %plot(rand(5));
end

% UIWAIT makes RRGMRES_DEMO wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = ...
    sym_rrgmres_demo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;
matrix_type=get(handles.matrix_type,'Value');
choice=get(handles.choice,'Value');
order=get(handles.order,'String'); order=str2double(order);
seed=get(handles.seed,'String'); seed=str2double(seed);
relerr=get(handles.relerror,'String'); relerr=str2double(relerr);
discrep_Value=get(handles.discrep_choice,'Value');
table_position = get(handles.errors,'Position');
total_column_width = (table_position(3) - 25);

if matrix_type==1
    set(handles.errors,...
        'ColumnName', ...
            {'Iterations', ...
            'Relative Residual Error', ...
            'Relative Solution Error'},...
        'ColumnWidth',...
            {'auto','auto','auto'}...
        )
    switch choice
        case 1
            [A,b,xexact]=baart(order);  %Baart Full
        case 2
            [A,b,xexact]=baart_alt(order);  %Baart Alt Full
        case 3
            [A,b,xexact]=shaw(order);   %Shaw
        case 4
            [A,b,xexact]=shaw_alt(order);   %Shaw Alt
        case 5
            [A,b,xexact]=deriv2(order,1);   %Deriv2 Linear
            xexact=xexact+1;
        case 6
            [A,b,xexact]=deriv2(order,2);   %Deriv2 Exp
            xexact=xexact+1;
        case 7
            [A,b,xexact]=deriv2(order,3);   %Deriv2 PW-Linear
            xexact=xexact+1;
        case 8
            [A,b,xexact]=deriv2_alt(order,1);   %Deriv2 Alt Linear
            xexact=xexact+1;
        case 9
            [A,b,xexact]=deriv2_alt(order,2);   %Deriv2 Alt Exp
            xexact=xexact+1;
        case 10
            [A,b,xexact]=phillips(order);   %Phillips
        case 11
            [A,b,xexact]=phillips_alt(order);   %Phillips Alt
    end
    b=A*xexact;
    randn('state',seed);
    err=randn(length(b),1); err=relerr*err*norm(b)/norm(err);
    b=b+err;
    
    if discrep_Value==1
        discrep=get(handles.discrep,'String'); discrep=str2double(discrep);
        if get(handles.discrep_value_choice,'Value') ~= 1
            delta = norm(err);
            discrep = delta*discrep;
        end
        
        [rrgmressol,~,~]=rrgmres_dp(A,b,discrep);
        for i = 1 : size(rrgmressol,2)
            %iteration
            rel_error(i,1) = i;
            %relative residual erro
            rel_error(i,2) = norm(A*rrgmressol(:,i)-b)/norm(b);
            %relative solution error
            rel_error(i,3) = norm(xexact-rrgmressol(:,i))/norm(xexact);
        end
        set(handles.errors,'Data',rel_error);
        plot(xexact,'k--'),hold,plot(rrgmressol(:,end),'k'),hold
        legend('exact','computed'); 
    else
        iterations=get(handles.spec_iter,'String'); 
        iterations=str2double(iterations);
        [rrgmressol]=rrgmres_iter(A,b,iterations);
        for i = 1 : size(rrgmressol,2)
            %iteration
            rel_error(i,1) = i;
            %relative residual error
            rel_error(i,2) = norm(A*rrgmressol(:,i)-b)/norm(b);
            %relative solution error
            rel_error(i,3) = norm(xexact-rrgmressol(:,i))/norm(xexact);
        end
        set(handles.errors,'Data',rel_error);
        plot(xexact,'k--'),hold,plot(rrgmressol(:,end),'k'),hold
        legend('exact','computed'); 
    end    
else
    %set(handles.errors,...
    %    'ColumnName', ...
    %        {'Iterations', ...
    %        'Relative Residual Error'}, ...
    %    'ColumnWidth',...
    %        {total_column_width/2, total_column_width/2}...
    %    )
    A=str2num(get(handles.user_A,'String'));
    b=str2num(get(handles.user_b,'String'));
    randn('state',seed);
    err=randn(length(b),1); err=relerr*err*norm(b)/norm(err);
    b=b+err;

    if discrep_Value==1
        discrep=get(handles.discrep,'String'); discrep=str2double(discrep);
        if get(handles.discrep_value_choice,'Value') ~= 1
            delta = norm(err);
            discrep = delta*discrep;
        end
        
        [rrgmressol,~,~]=rrgmres_dp(A,b,discrep);
        
        for i = 1 : size(rrgmressol,2)
            %iteration
            rel_error(i,1) = i; 
            %relative residual erorr
            rel_error(i,2) = norm(A*rrgmressol(:,i)-b)/norm(b);
        end
        set(handles.errors,'Data',rel_error);

        plot(rrgmressol(:,end),'k')
        legend('computed'); 
    else
        iterations=get(handles.spec_iter,'String'); 
        iterations=str2double(iterations);
        [rrgmressol]=rrgmres_iter(A,b,iterations);
        for i = 1 : size(rrgmressol,2)
            %iteration
            rel_error(i,1) = i;
            %relative residual erro
            rel_error(i,2) = norm(A*rrgmressol(:,i)-b)/norm(b);
        end
        set(handles.errors,'Data',rel_error);
    
        plot(rrgmressol(:,end),'k')
        legend('computed'); 
    end    
end



end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end
end
% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)
end
% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

end
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') 
%        returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from 
%        popupmenu1

end
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', ...
    'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
end


function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of order as text
%        str2double(get(hObject,'String')) returns contents of order as a 
%        double

end
% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in choice.
function choice_Callback(hObject, eventdata, handles)
% hObject    handle to choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) 
%       returns choice contents as cell array 
%       contents{get(hObject,'Value')} returns selected item from choice

end
% --- Executes during object creation, after setting all properties.
function choice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function seed_Callback(hObject, eventdata, handles)
% hObject    handle to seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seed as text
%        str2double(get(hObject,'String')) returns contents of seed as a 
%        double

end
% --- Executes during object creation, after setting all properties.
function seed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function relerror_Callback(hObject, eventdata, handles)
% hObject    handle to relerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of relerror as text
%        str2double(get(hObject,'String')) returns contents of relerror as 
%        a double

end
% --- Executes during object creation, after setting all properties.
function relerror_CreateFcn(hObject, eventdata, handles)
% hObject    handle to relerror (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function discrep_Callback(hObject, eventdata, handles)
% hObject    handle to discrep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of discrep as text
%        str2double(get(hObject,'String')) returns contents of discrep as 
%        a double
end

% --- Executes during object creation, after setting all properties.
function discrep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to discrep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
% --------------------------------------------------------------------
function uipanel1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


end
function spec_iter_Callback(hObject, eventdata, handles)
% hObject    handle to spec_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spec_iter as text
%        str2double(get(hObject,'String')) returns contents of spec_iter 
%        as a double


% --- Executes during object creation, after setting all properties.
end
function spec_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spec_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end
% --- Executes during object creation, after setting all properties.
function report_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to report_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
% --- Executes during object creation, after setting all properties.
function report_relres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to report_relres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
% --- Executes during object creation, after setting all properties.
function report_relsol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to report_relsol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end



function user_A_Callback(hObject, eventdata, handles)
% hObject    handle to user_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of user_A as text
%        str2double(get(hObject,'String')) returns contents of user_A as a 
%        double

end
% --- Executes during object creation, after setting all properties.
function user_A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function user_b_Callback(hObject, eventdata, handles)
% hObject    handle to user_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of user_b as text
%        str2double(get(hObject,'String')) returns contents of user_b as a 
%        double
end

% --- Executes during object creation, after setting all properties.
function user_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to user_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function report_iter_Callback(hObject, eventdata, handles)
% hObject    handle to report_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of report_iter as text
%        str2double(get(hObject,'String')) returns contents of report_iter 
%        as a double
end


function report_relres_Callback(hObject, eventdata, handles)
% hObject    handle to report_relres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of report_relres as text
%        str2double(get(hObject,'String')) returns contents of 
%        report_relres as a double

end

function report_relsol_Callback(hObject, eventdata, handles)
% hObject    handle to report_relsol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of report_relsol as text
%        str2double(get(hObject,'String')) returns contents of 
%        report_relsol as a double
end
