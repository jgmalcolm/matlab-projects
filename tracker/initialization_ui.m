function varargout = initialization_ui(varargin)
% INITIALIZATION_UI M-file for initialization_ui.fig
%      INITIALIZATION_UI, by itself, creates a new INITIALIZATION_UI or raises the existing
%      singleton*.
%
%      H = INITIALIZATION_UI returns the handle to a new INITIALIZATION_UI or the handle to
%      the existing singleton*.
%
%      INITIALIZATION_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INITIALIZATION_UI.M with the given input arguments.
%
%      INITIALIZATION_UI('Property','Value',...) creates a new INITIALIZATION_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before initialization_ui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to initialization_ui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help initialization_ui

% Last Modified by GUIDE v2.5 22-Jan-2006 18:00:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @initialization_ui_OpeningFcn, ...
                   'gui_OutputFcn',  @initialization_ui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before initialization_ui is made visible.
function initialization_ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to initialization_ui (see VARARGIN)

for i=1:2:length(varargin)
    if( strcmp( varargin{i}, 'first_image' ) == 1 )
        handles.first_image = varargin{2};   
        imagesc( handles.first_image );  axis image;   
        colormap gray;
        set(get(handles.axes1,'Children'),'ButtonDownFcn', ...
            'initialization_ui(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
    end
end
handles.class_selection_state = 0;
handles.UseDynamics = 0;
handles.class_points = [];
% Choose default command line output for initialization_ui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes initialization_ui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = initialization_ui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in AddClassButton.
function AddClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if( handles.class_selection_state == 0 )
    handles.class_selection_state = 1; 
    handles.current_points = [];
    set(handles.AddClassButton,'Enable','off');
end

% Choose default command line output for initialization_ui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in DoneWithClassButton.
function DoneWithClassButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneWithClassButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( handles.class_selection_state == 1 ) 
  if( isfield( handles, 'class_points' ) == 0 ) 
    handles.class_points{1} = handles.current_points;
  else
    l = length( handles.class_points );
    handles.class_points{l+1} = handles.current_points; 
  end
  handles.current_points = [];
  handles.class_selection_state = 0;
  set(handles.ClassCounter,'String',num2str(length(handles.class_points)));
  if( length(handles.class_points) == 1 )
    set(handles.DoneButton,'Enable','on');
  end
  set(handles.DoneWithClassButton,'Enable','off');
  set(handles.AddClassButton,'Enable','on');
end

% Update handles structure
guidata(hObject, handles);

%function that sets the use dynamics flag
function UseDynamicsButton_Callback(hObject, eventdata, handles)
 handles.UseDynamics = 1;
 set(handles.UseDynamicsButton,'Enable','off');
 
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global initialization_ui_output;

if(isempty(handles.class_points))
    delete(gcf);
    return;
end
initialization_ui_output.class_points = handles.class_points;
initialization_ui_output.UseDynamics = handles.UseDynamics;
%handles.class_points
delete(gcf);


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( handles.class_selection_state == 1 ) 
  the_point = ...
      get(handles.axes1,'CurrentPoint');
  handles.current_points = [ handles.current_points ; ...
                      [ the_point(1,1),the_point(1,2) ] ];
  num_points = size( handles.current_points, 1 );
  hold on; 
%     plot( handles.current_points(end,1), ...
% 	  handles.current_points(end,2), 'ro' );
  if( num_points > 1 ) 
    plot( handles.current_points((num_points-1):num_points,1), ...
	  handles.current_points((num_points-1):num_points,2), 'r-' );
  end
  hold off;
  if( strcmp( get(handles.DoneWithClassButton,'Enable'), 'off' ) == 1 )
    set(handles.DoneWithClassButton,'Enable','on');
  end
end

% Update handles structure
guidata(hObject, handles);
