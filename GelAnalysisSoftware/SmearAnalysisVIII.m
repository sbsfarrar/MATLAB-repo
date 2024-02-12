function varargout = SmearAnalysisVIII(varargin)
% SMEARANALYSISVIII M-file for SmearAnalysisVIII.fig
%      SMEARANALYSISVIII, by itself, creates a new SMEARANALYSISVIII or raises the existing
%      singleton*.
%
%      H = SMEARANALYSISVIII returns the handle to a new SMEARANALYSISVIII or the handle to
%      the existing singleton*.
%
%      SMEARANALYSISVIII('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMEARANALYSISVIII.M with the given input arguments.
%
%      SMEARANALYSISVIII('Property','Value',...) creates a new SMEARANALYSISVIII or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SmearAnalysisVIII_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SmearAnalysisVIII_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SmearAnalysisVIII

% Last Modified by GUIDE v2.5 27-May-2016 11:14:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SmearAnalysisVIII_OpeningFcn, ...
    'gui_OutputFcn',  @SmearAnalysisVIII_OutputFcn, ...
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


% --- Executes just before SmearAnalysisVIII is made visible.
function SmearAnalysisVIII_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SmearAnalysisVIII (see VARARGIN)

% Choose default command line output for SmearAnalysisVIII
handles.output = hObject;
handles.x=[];
handles.smear=[];
handles.smear_filt=[];
handles.Gel_filt=[];
handles.expfit={};
handles.IratioSeries=[];
handles.Bandfit={};
handles.Areas=[];
handles.pnum=0;
handles.PeakExport={};
handles.Iall={};
handles.Xall={};
handles.LeftBG=[];
handles.RightBG=[];
% Update handles structure
handles.BfitOption=1
handles.popup=1;
handles.boxdraw=0;
handles.AreasAll=[];
handles.AreasNormAll=[];
guidata(hObject, handles);

% UIWAIT makes SmearAnalysisVIII wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SmearAnalysisVIII_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadImage.
function LoadImage_Callback(hObject, ~, handles)
% hObject    handle to LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.image);
[filename,pathname] = uigetfile('*.*','Select file');
f=filesep;
handles.Gel=imread([pathname,f,filename]);

handles.tfile=Tiff([pathname,f,filename],'r');
%axes(handles.image);
imshow(handles.Gel);
guidata(hObject, handles);

% --- Executes on button press in ROI.
function ROI_Callback(hObject, eventdata, handles)
% hObject    handle to ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.image);
%handles.box=impoly('Closed',true);
%position = wait(handles.box);


handles.box=imrect;
handles.coord=round(getPosition(handles.box));

%Make sure all the axes line up


guidata(hObject, handles);

% --- Executes on button press in PlotROI.
function PlotROI_Callback(hObject, eventdata, handles)
% hObject    handle to PlotROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.coord=round(getPosition(handles.box));


%Assigne Left, Right, Top and Bottom
T=handles.coord(2);
B=handles.coord(2)+handles.coord(4);
L=handles.coord(1);
R=handles.coord(1)+handles.coord(3);
% %if get(handles.MidEnable,'Value')==1
%     subImg1=handles.Gel(T:B,L:midp)
%     subImg2=handles.Gel(T:B,midp:R);
%     subImg1=subImg1';
%     subImg2=subImg2';
%     handles.smear=sum(subImg2)-sum(subImg1);
% else
%Thresh=str2num(get(handles.Thresh_Median,'string'));

subImg=handles.Gel(T:B,L:R);
subImg=subImg;

SpotRemove=[];
handles.smear=[];
handles.smear_filtMode=[];
handles.smear_filtMedian=[];
handles.smear_filtMean_of_Subpopulation=[];
handles.smear_filt=[];
handles.smear_sum=[];
%     subImg=subImg';
for i=1:length(subImg(:,1))
    Ir=subImg(i,:);
    Ir=cast(Ir,'double');
    IrS=Ir;
    IrMSP=Ir;
    handles.smear_filt(i)=median(Ir);
    %         f=find(Ir<Thresh*median(Ir));
    %         IrS(f)=[];
    %         NewThresh=3*std(cast(IrS,'double'));
    %         f=find(Ir<median(Ir)-NewThresh);
    %         Ir(f)=[];
    handles.smear(i)=mean(Ir);
    handles.smear_sum(i)=sum(Ir);
    handles.smear_filtMode(i)=mode(round(Ir/100))*100;
    handles.smear_filtMedian(i)=median(Ir);
    
    sortedIr = sort(IrMSP,'ascend');
    len = length(IrMSP);
    Percentile_range_start = 0.25;
    Percentile_range_end = 0.75;
    handles.smear_filtMean_of_Subpopulation(i) = mean(sortedIr(round(Percentile_range_start*len):round(Percentile_range_end*len)));
    
    
end
%     handles.smear=sum(subImg);
%     if isempty(handles.Gel_filt)==0
%         subImg_filt=handles.Gel_filt(T:B,L:R);
%         subImg_filt=subImg_filt';
%         handles.smear_filt=sum(subImg_filt);
%         handles.subImg_filt=subImg_filt';
%         handles.smear_filt=max(handles.smear)-handles.smear_filt;
%
%
%     end

% end
%handles.subImg=subImg';
handles.smear_filtMedian=max(handles.smear_filtMedian)-handles.smear_filtMedian;
handles.smear_filtMode=max(handles.smear_filtMode)-handles.smear_filtMode;
handles.smear=max(handles.smear)-handles.smear;
handles.smear_filt=max(handles.smear_filt)-handles.smear_filt;
handles.smear_sum=max(handles.smear_sum)-handles.smear_sum;
handles.smear_filtMean_of_Subpopulation = max(handles.smear_filtMean_of_Subpopulation) - handles.smear_filtMean_of_Subpopulation;
check=get(handles.MedianFilter,'Value');
if check==0
    handles.smear_filt=handles.smear;
end
check=get(handles.SumIntensity,'Value');

if check==1
    handles.smear_filt=handles.smear_sum;
    handles.smear=handles.smear_sum;
end
% max(handles.smear);
% handles.smear=handles.smear./max(handles.smear);
L=length(handles.smear);
axes(handles.Rplot);
handles.x=[1:L];
% plot(handles.x,handles.smear,'b',handles.x,handles.smear_filtMedian,'r',handles.x,handles.smear_filtMode,'g',handles.x,handles.smear_filtMean_of_Subpopulation,'bo');
plot(handles.x,handles.smear_filt,'b',handles.x,handles.smear,'r');
axis tight;


% if handles.pnum==0
%     handles.c1=MakeCursor(handles.Rplot,handles.x(round(L/3)),handles.smear(round(L/3)),'Red');
%     handles.c2=MakeCursor(handles.Rplot,handles.x(round(2*L/3)),handles.smear(round(2*L/3)),'Green');
%     handles.pnum=handles.pnum+1;
% else
%     handles.pnum=handles.pnum+1;
% end

if ~isempty(handles.BGcoord)
    
   handles.bgbox=impoly(gca,handles.BGcoord,'Closed',false); 
end
guidata(hObject, handles);




function maxX_Callback(hObject, eventdata, handles)
% hObject    handle to maxX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxX as text
%        str2double(get(hObject,'String')) returns contents of maxX as a double


% --- Executes during object creation, after setting all properties.
function maxX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FitGel.
function FitGel_Callback(hObject, eventdata, handles)
% hObject    handle to FitGel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



guidata(hObject, handles);
p1=handles.c1.position;
xmin=round(p1(1));
p2=handles.c2.position;
xmax=round(p2(1));
% xmxstr=get(handles.maxX,'String');
% xmnstr=get(handles.minX,'String');
dvstr=get(handles.dv,'String');
dv=str2num(dvstr);
% xmax=str2num(xmxstr);
% xmin=str2num(xmnstr);
xfit=xmin:xmax;
smfit=handles.smear(xmin:xmax);
smstart=smfit(1);
smfit=smfit;
tmax=xmax/dv;
sigmastr=get(handles.sigma,'String');
sigma=str2num(sigmastr);
xstart=xfit(1);
xfit=xfit-xstart;
Aguess=max(handles.smear);
switch handles.popup
    
    case 1
        %Standard Exponential Fit
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0,],...
            'Upper',[Aguess*2,3*xmax,],...
            'Startpoint',[Aguess xmax/2 ]);
        f = fittype('A*exp(-x/tx)','options',s);
        [fitdata,fitstats]=fit(xfit',smfit(end:-1:1)',f);
        yp=fitdata.A*exp(-(xfit)./(fitdata.tx));
        handles.sfit=fitdata;
    case 2
        %Exponetial Convolution Fit
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0, -30],...
            'Upper',[2*Aguess,3*xmax, xmax],...
            'Startpoint',[Aguess xmax xmax/2 ]);
        
        f = fittype('A*exp(-(x-xs)/tx)*(erf(1/2*(x-xs)*2^(1/2)/s)+1)','options',s,'problem','s');
        [fitdata,fitstats]=fit(xfit',smfit(end:-1:1)',f,'problem',handles.bfit.s);
        handles.sfit=fitdata;
        yp=fitdata.A*exp(-(xfit-fitdata.xs)./(fitdata.tx)).*(erf(1/2*(xfit-fitdata.xs)*2^(1/2)/fitdata.s)+1);
    case 3
        %Exponential Convolution + Gaussian
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0,-30 0 0 xmax/2 ],...
            'Upper',[2*Aguess,3*xmax, xmax 3*handles.bfit.s 3*max(smfit) xmax ],...
            'Startpoint',[Aguess xmax xmax/2 handles.bfit.s (max(smfit)-smfit(end)) handles.bfit.xo ]);
        
        f = fittype('A*exp(-(x-xs)/tx)*(erf(1/2*(x-xs)*2^(1/2)/s)+1)+B*exp(-(x-xm)^2/2/s^2)','options',s);
        [fitdata,fitstats]=fit(xfit',smfit(end:-1:1)',f);
        yp=fitdata.A*exp(-(xfit-fitdata.xs)./(fitdata.tx)).*(erf(1/2*(xfit-fitdata.xs)*2^(1/2)/fitdata.s)+1)+smstart+fitdata.B*exp(-(xfit-fitdata.xm).^2/2/fitdata.s^2);
end





axes(handles.Rplot);
hold on;

%yp=fitdata.A*exp(-(xfit)./(fitdata.tx)).*(erf(1/2*(xfit)*2^(1/2)/sigma)+1)+smstart;
dv=str2double(get(handles.dv,'String'));
set(handles.tau,'String',num2str(fitdata.tx/dv));

plot(xfit+xstart,yp(end:-1:1),'r');
%fitdata
hold off;
guidata(hObject, handles);
f1=figure;
plot(handles.x,handles.smear);
hold on;plot(xfit+xstart,yp(end:-1:1),'r');

handles.sfit=fitdata;


function minX_Callback(hObject, eventdata, handles)
% hObject    handle to minX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minX as text
%        str2double(get(hObject,'String')) returns contents of minX as a double


% --- Executes during object creation, after setting all properties.
function minX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dv_Callback(hObject, eventdata, handles)
% hObject    handle to dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dv as text
%        str2double(get(hObject,'String')) returns contents of dv as a double


% --- Executes during object creation, after setting all properties.
function dv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function midp_Callback(hObject, eventdata, handles)
% hObject    handle to midp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of midp as text
%        str2double(get(hObject,'String')) returns contents of midp as a double


% --- Executes during object creation, after setting all properties.
function midp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to midp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MidEnable.
function MidEnable_Callback(hObject, eventdata, handles)
% hObject    handle to MidEnable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MidEnable



function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma as text
%        str2double(get(hObject,'String')) returns contents of sigma as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bstart_Callback(hObject, eventdata, handles)
% hObject    handle to Bstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bstart as text
%        str2double(get(hObject,'String')) returns contents of Bstart as a double


% --- Executes during object creation, after setting all properties.
function Bstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bend_Callback(hObject, eventdata, handles)
% hObject    handle to Bend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bend as text
%        str2double(get(hObject,'String')) returns contents of Bend as a double


% --- Executes during object creation, after setting all properties.
function Bend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FitBand.
function FitBand_Callback(hObject, eventdata, handles)
% hObject    handle to FitBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);


p1=handles.c1.position;
xmin=round(p1(1));
p2=handles.c2.position;
xmax=round(p2(1));

% xmxstr=get(handles.Bend,'String');
% xmnstr=get(handles.Bstart,'String');
% xmax=str2num(xmxstr);
% xmin=str2num(xmnstr);
xfit=xmin:xmax;
smfit=handles.smear(xmin:xmax);


switch handles.BfitOption
    
    case 1
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0,0],...
            'Upper',[2*max(handles.smear),xmax,4*xmax],...
            'Startpoint',[max(handles.smear) mean(xfit) (xmax-xmin)/2 ] );
        f = fittype('A*exp(-(x-xo)^2/2/s^2)','options',s);
        [handles.bfit,handles.bfits]=fit(xfit',smfit',f);
        axes(handles.Rplot);
        xfit=xmin-10:0.1:xmax+4;
        set(handles.sigma,'String',num2str(handles.bfit.s));
        yp=handles.bfit.A*exp(-(xfit-handles.bfit.xo).^2/2/handles.bfit.s^2);
        
        
    case 2
        s=str2double(get(handles.sigma,'String'));
        b1=str2double(get(handles.B1Guess,'String'));
        b2=str2double(get(handles.B2Guess,'String'));
        db=str2double(get(handles.dpMax,'String'));
        f = fittype('gauss2');
        options=fitoptions('gauss2');
        options.Lower=[0 b1-db 0 0 b2-db 0];
        options.Upper=[1.1*max(smfit) b1+db 1.2*s 1.1*max(smfit) b2+db 1.2*s];
        [bfit,handles.bfits]=fit(xfit',smfit',f);
        xfit=xmin-10:0.1:xmax+4;
        yp=bfit.a1*exp(-((xfit-bfit.b1)/bfit.c1).^2) + bfit.a2*exp(-((xfit-bfit.b2)/bfit.c2).^2);
        handles.bfit=bfit;
        
end
hold on;
plot(xfit,yp,'r');
hold off;
guidata(hObject, handles);
%         f1=figure;
%         plot(handles.smear);
%         hold on;
%         plot(xfit,yp,'r');
%         hold off;
%
guidata(hObject, handles);

function minS_Callback(hObject, eventdata, handles)
% hObject    handle to minS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minS as text
%        str2double(get(hObject,'String')) returns contents of minS as a double


% --- Executes during object creation, after setting all properties.
function minS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxS_Callback(hObject, eventdata, handles)
% hObject    handle to maxS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxS as text
%        str2double(get(hObject,'String')) returns contents of maxS as a double


% --- Executes during object creation, after setting all properties.
function maxS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SubtractBand.
function SubtractBand_Callback(hObject, eventdata, handles)
% hObject    handle to SubtractBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
p1=handles.c1.position;
xmin=round(p1(1));
p2=handles.c2.position;
xmax=round(p2(1));
% xmxstr=get(handles.maxS,'String');
% xmnstr=get(handles.minS,'String');
% xmax=str2num(xmxstr);
% xmin=str2num(xmnstr);
xS=xmin:xmax;
yS=handles.bfit.A*exp(-(xS-handles.bfit.xo).^2/2/handles.bfit.s^2);
handles.smear(xmin:xmax)=handles.smear(xmin:xmax)-yS;
axes(handles.Rplot);
plot(handles.smear);
guidata(hObject, handles);


% --- Executes on button press in LoadBandFit.
function LoadBandFit_Callback(hObject, eventdata, handles)
% hObject    handle to LoadBandFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

handles.expfit{end+1}=handles.bfit;
guidata(hObject, handles);
% --- Executes on button press in LoadSmearFit.
function LoadSmearFit_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSmearFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

handles.expfit{end+1}=handles.sfit;
guidata(hObject, handles);
% --- Executes on button press in ClearFitExp.
function ClearFitExp_Callback(hObject, eventdata, handles)
% hObject    handle to ClearFitExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.expfit={};
guidata(hObject, handles);
% --- Executes on button press in ExpFit.
function ExpFit_Callback(hObject, eventdata, handles)
% hObject    handle to ExpFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
assignin('base', 'TempExpFitData', handles.expfit);
evalin('base','ExpFitData{end+1}=TempExpFitData;');

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FitPopUp.
function FitPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to FitPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.popup=get(hObject,'Value');
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns FitPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FitPopUp


% --- Executes during object creation, after setting all properties.
function FitPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FitPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function IntensitySliderMax_Callback(hObject, eventdata, handles)
% hObject    handle to IntensitySliderMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% ImgTemp=handles.Gel;
% f=find(ImgTemp>=Mx);
% ImgTemp(f)=Mx;
handles.Mx=get(handles.IntensitySliderMax,'Value');

axes(handles.image);
imshow(handles.Gel,[handles.Mn handles.Mx]);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function IntensitySliderMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntensitySliderMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',3);
set(hObject,'Min',1);
set(hObject,'Max',2^16);
set(hObject,'Value',2^16);
handles.Mx=2^16;
guidata(hObject, handles);
% --- Executes on slider movement.
function IntensitySliderMin_Callback(hObject, eventdata, handles)
% hObject    handle to IntensitySliderMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.Mn=get(handles.IntensitySliderMin,'Value');

axes(handles.image);
imshow(handles.Gel,[handles.Mn handles.Mx]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function IntensitySliderMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IntensitySliderMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',3);
set(hObject,'Min',1);
set(hObject,'Max',2^16);
handles.Mn=3;
guidata(hObject, handles);



function dist_Callback(hObject, eventdata, handles)
% hObject    handle to dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist as text
%        str2double(get(hObject,'String')) returns contents of dist as a double


% --- Executes during object creation, after setting all properties.
function dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Ratio.
function Ratio_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d=str2num(get(handles.dist,'String'));
index=round(handles.bfit.xo-d);
handles.Iratio=handles.smear(index)/(handles.bfit.A+handles.smear(index));
guidata(hObject, handles);


% --- Executes on button press in LoadRatio.
function LoadRatio_Callback(hObject, eventdata, handles)
% hObject    handle to LoadRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.IratioSeries(end+1)=handles.Iratio;
guidata(hObject, handles);
% --- Executes on button press in ExpRatio.
function ExpRatio_Callback(hObject, eventdata, handles)
% hObject    handle to ExpRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
assignin('base', 'IratioSeries', handles.IratioSeries)

% --- Executes on button press in ClearRatio.
function ClearRatio_Callback(hObject, eventdata, handles)
% hObject    handle to ClearRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.IratioSeries=[];
guidata(hObject, handles);


% --- Executes on button press in DrawBG.
function DrawBG_Callback(hObject, eventdata, handles)
% hObject    handle to DrawBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guidata(hObject, handles);
Nmax=str2double(get(handles.vnum,'String'));
% x1=1;
L=length(handles.smear);
% x2=L/3;
% x3=2*L/3;
% x4=L;
Yall=max(handles.smear)/2;
coords(1,1:2)=[1,Yall];
d=Nmax-1;
for i=2:Nmax-1
    coords=[coords;[L*(i-1)/d,Yall]];
end
coords=[coords;[L,Yall]];


%[x1,Yall; x2,Yall; x3,Yall; x4,Yall]
axes(handles.Rplot);
if handles.boxdraw>0
    handles.bgbox.delete;
end
handles.boxdraw=handles.boxdraw+1;
handles.bgbox=impoly(gca,[coords],'Closed',false);
guidata(hObject, handles);

% --- Executes on button press in SubtractBG.
function SubtractBG_Callback(hObject, eventdata, handles)
% hObject    handle to SubtractBG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);


state=get(handles.BGTemplate,'Value');
if state==0
    handles.BGcoord=round(getPosition(handles.bgbox));
    
end
BGcoord=handles.BGcoord;

BGcoord(:,2)=handles.smear(BGcoord(:,1));
Nmax=length(BGcoord(:,1));
Ibgr=[];

for i=2:Nmax
    
    %     x2=BGcoord(i,1);y2=BGcoord(i,2);
    %     x1=BGcoord(i-1,1);y1=BGcoord(i-1,2);
    %
    x2=BGcoord(i,1);y2=handles.smear_filt(x2);
    x1=BGcoord(i-1,1);y1=handles.smear_filt(x1);
    
    m=(y2-y1)/(x2-x1);
    b=y2-m*x2;
    xr=x1:x2;
    IbgrT=xr.*m+b;
    if i>2
        x1=1+x1;
        Ibgr=[Ibgr,IbgrT(2:end)];
    else
        Ibgr=IbgrT;
    end
    
    
end
handles.smear_filt(BGcoord(1,1):BGcoord(end,1))=handles.smear_filt(BGcoord(1,1):BGcoord(end,1))-Ibgr;
plot(handles.x,handles.smear_filt);
handles.bgbox=impoly(gca,handles.BGcoord,'Closed',false);
L=length(handles.smear);
% handles.c1=MakeCursor(handles.Rplot,handles.x(round(L/3)),handles.smear(round(L/3)),'Red');
% handles.c2=MakeCursor(handles.Rplot,handles.x(round(2*L/3)),handles.smear(round(2*L/3)),'Green');
BGcorrect.I=Ibgr;
BGcorrect.x1=BGcoord(1,1);
BGcorrect.x2=BGcoord(1,end);
handles.BGcorrect=BGcorrect;

AreaCheck=get(handles.AreaCheck,'Value');

if AreaCheck==1
    
    handles.Areas=AreaFromPoint(BGcoord,handles.smear_filt);
    Total=sum(handles.Areas);
    handles.AreasNorm=handles.Areas./Total;
end

guidata(hObject, handles);



function pthresh_Callback(hObject, eventdata, handles)
% hObject    handle to pthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pthresh as text
%        str2double(get(hObject,'String')) returns contents of pthresh as a double


% --- Executes during object creation, after setting all properties.
function pthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindPeaks.
function FindPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to FindPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

minp=str2double(get(handles.pthresh,'String'))/100;
mins=str2double(get(handles.minW,'String'))/100;
[peaks, PFWHH,PExt]=mspeaks(handles.x,handles.smear,'HeightFilter',max(handles.smear)*minp,'ShowPlot',false);
%wmax=max(PExt(:,2)-PExt(:,1));
%[peaks]=mspeaks(handles.x,handles.smear,'HeightFilter',max(handles.smear)*minp,'FWHHFilter', mins*wmax);

axes(handles.Rplot);
hold off;
plot(handles.smear);hold on;
plot(peaks(:,1),peaks(:,2),'go');hold off;
handles.peaks=peaks;
guidata(hObject, handles);



function minW_Callback(hObject, eventdata, handles)
% hObject    handle to minW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minW as text
%        str2double(get(hObject,'String')) returns contents of minW as a double


% --- Executes during object creation, after setting all properties.
function minW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vnum_Callback(hObject, eventdata, handles)
% hObject    handle to vnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vnum as text
%        str2double(get(hObject,'String')) returns contents of vnum as a double


% --- Executes during object creation, after setting all properties.
function vnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FitPeaks.
function FitPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to FitPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
s=str2double(get(handles.sigma,'String'));
dmin=str2double(get(handles.MinD,'String'));
range=str2double(get(handles.MaxR,'String'));
FitInfo=GetFitInfo(handles.peaks,dmin);
PeakFitData=FitPeaks(FitInfo,handles.smear,range);
axes(handles.Rplot);
hold off;
plot(handles.smear);
hold on;
for i=1:length(PeakFitData)
    plot(PeakFitData{i});
end
hold off;

L=length(handles.smear);
% handles.c1=MakeCursor(handles.Rplot,handles.x(round(L/3)),handles.smear(round(L/3)),'Red');
% handles.c2=MakeCursor(handles.Rplot,handles.x(round(2*L/3)),handles.smear(round(2*L/3)),'Green');


IntensitySum=SumIntensityFit(PeakFitData,FitInfo);
handles.PeakFitData=PeakFitData;
handles.IntensitySum=IntensitySum;
handles.FitInfo=FitInfo;
guidata(hObject, handles);
function MinD_Callback(hObject, eventdata, handles)
% hObject    handle to MinD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinD as text
%        str2double(get(hObject,'String')) returns contents of MinD as a double


% --- Executes during object creation, after setting all properties.
function MinD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxR_Callback(hObject, eventdata, handles)
% hObject    handle to MaxR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxR as text
%        str2double(get(hObject,'String')) returns contents of MaxR as a double


% --- Executes during object creation, after setting all properties.
function MaxR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadPeaks.
function LoadPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
dims=size(handles.Areas);
colA=dims(2);
colIS=length(handles.IntensitySum);
if colA>=colIS
    istart=colA-colIS+1;
    handles.Areas(end+1,istart:end)=handles.IntensitySum;
else
    rowA=dims(1);
    cadd=colIS-colA;
    z=zeros(rowA,cadd);
    handles.Areas=[z,handles.Areas];
    handles.Areas=[handles.Areas;handles.IntensitySum];
end



Bandfit.PeakFitData=handles.PeakFitData;
Bandfit.FitInfo=handles.FitInfo;
Bandfit.BGcorrect=handles.BGcorrect;
Bandfit.Iplot=handles.smear;
Bandfit.subImg=handles.subImg;
handles.Bandfit{end+1}=Bandfit;
guidata(hObject, handles);
% --- Executes on button press in ExpPeaks.
function ExpPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to ExpPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
assignin('base', 'BandFit', handles.Bandfit);
assignin('base', 'Areas', handles.Areas);
% --- Executes on button press in ClearPeaks.
function ClearPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to ClearPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.Areas=[];
handles.Bandfit={};
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function Rplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Rplot



function GelTime_Callback(hObject, eventdata, handles)
% hObject    handle to GelTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GelTime as text
%        str2double(get(hObject,'String')) returns contents of GelTime as a double


% --- Executes during object creation, after setting all properties.
function GelTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GelTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GetdV.
function GetdV_Callback(hObject, eventdata, handles)
% hObject    handle to GetdV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
p1=handles.c1.position;
xmin=round(p1(1));
p2=handles.c2.position;
xmax=round(p2(1));
dx=abs(xmax-xmin);
t=str2num(get(handles.GelTime,'String'));
handles.dvN=dx/t;
set(handles.dv,'String',num2str(handles.dvN));


% --- Executes on button press in AddProbes.
function AddProbes_Callback(hObject, eventdata, handles)
% hObject    handle to AddProbes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
L=length(handles.smear);
handles.cursorMode=datacursormode;
handles.c1=MakeCursor(handles.Rplot,handles.x(round(L/3)),handles.smear(round(L/3)),'Red',handles.cursorMode);
handles.c2=MakeCursor(handles.Rplot,handles.x(round(2*L/3)),handles.smear(round(2*L/3)),'Green',handles.cursorMode);
guidata(hObject, handles);

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CropImage.
function CropImage_Callback(hObject, eventdata, handles)
% hObject    handle to CropImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.coord=round(getPosition(handles.box));
%midPs=get(handles.midp,'String');
%midp=str2num(midPs);
xmin=handles.coord(1);
ymin=handles.coord(2);
w=handles.coord(3);
h=handles.coord(4);
%Make sure all the axes line up
% handles.coord(end,1)=handles.coord(1,1);
% handles.coord(2,1)=handles.coord(3,1);
% handles.coord(1,2)=handles.coord(2,2);
% handles.coord(3,2)=handles.coord(4,2);
% setPosition(handles.box,handles.coord);
%Assigne Left, Right, Top and Bottom
L=xmin;
R=xmin+w;
T=ymin;
B=ymin+h;
% %if get(handles.MidEnable,'Value')==1
%     subImg1=handles.Gel(T:B,L:midp)
%     subImg2=handles.Gel(T:B,midp:R);
%     subImg1=subImg1';
%     subImg2=subImg2';
%     handles.smear=sum(subImg2)-sum(subImg1);
% else
subImg=handles.Gel(T:B,L:R);
handles.Gel=subImg;
axes(handles.image);
imshow(handles.Gel);
handles.box.delete;
guidata(hObject, handles);


% --- Executes on button press in RemoveNoise.
function RemoveNoise_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
p1=handles.c1.position;
x1=round(p1(1));y1=p1(2);
p2=handles.c2.position;
x2=round(p2(1));y2=p2(2);


% m=(y2-y1)/(x2-x1);
% b=y2-m*x2;
% x=x1:x2;
% y=m*x+b;
handles.smear_filt(x1:x2)=[];
handles.x(x1:x2)=[];
axes(handles.Rplot);
plot(handles.x,handles.smear_filt);
guidata(hObject, handles);
L=length(handles.smea_filt);

guidata(hObject, handles);


% --- Executes on button press in InvertImage.
function InvertImage_Callback(hObject, eventdata, handles)
% hObject    handle to InvertImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.Gel=imcomplement(handles.Gel);
axes(handles.image);
imshow(handles.Gel);
guidata(hObject, handles);


% --- Executes on button press in LockH.
function LockH_Callback(hObject, eventdata, handles)
% hObject    handle to LockH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LockH
guidata(hObject, handles);
state=get(handles.LockH,'Value');
s=size(handles.Gel);
coord=getPosition(handles.box);
switch(state)
    
    case(1)
        fcn=makeConstrainToRectFcn('imrect',[1 s(2)], [coord(2)+coord(4) coord(2)+coord(4)]);
        setPositionConstraintFcn(handles.box,fcn);
        
    case(0)
        fcn=makeConstrainToRectFcn('imrect',[1 s(2)], [1 s(1)]);
        setPositionConstraintFcn(handles.box,fcn);
        
end



% --------------------------------------------------------------------
function uitoggletool3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in BandOption.
function BandOption_Callback(hObject, eventdata, handles)
% hObject    handle to BandOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BandOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BandOption

guidata(hObject, handles);
handles.BfitOption=get(hObject,'Value');
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function BandOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BandOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B2Guess_Callback(hObject, eventdata, handles)
% hObject    handle to B2Guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B2Guess as text
%        str2double(get(hObject,'String')) returns contents of B2Guess as a double


% --- Executes during object creation, after setting all properties.
function B2Guess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B2Guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B1Guess_Callback(hObject, eventdata, handles)
% hObject    handle to B1Guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B1Guess as text
%        str2double(get(hObject,'String')) returns contents of B1Guess as a double


% --- Executes during object creation, after setting all properties.
function B1Guess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B1Guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dpMax_Callback(hObject, eventdata, handles)
% hObject    handle to dpMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dpMax as text
%        str2double(get(hObject,'String')) returns contents of dpMax as a double


% --- Executes during object creation, after setting all properties.
function dpMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dpMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function tau_Callback(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau as text
%        str2double(get(hObject,'String')) returns contents of tau as a double


% --- Executes during object creation, after setting all properties.
function tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in PlotSeperate.
function PlotSeperate_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSL>
% hObject    handle to PlotSeperate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f1=figure
figure(f1);
plot(handles.x,handles.smear)


% --- Executes on button press in ExpIntensity.
function ExpIntensity_Callback(hObject, eventdata, handles) %#ok<INUSD>
% hObject    handle to ExpIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);

assignin('base', 'TempXall', handles.Xall);
assignin('base', 'TempIall', handles.Iall);
assignin('base','AreasAll',handles.AreasAll);
assignin('base','AreasNormAll',handles.AreasNormAll);

% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
handles.PeakExport{end+1}=handles.peaks;
guidata(hObject, handles);

% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
assignin('base', 'PeakData', handles.PeakExport);

% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
handles.PeakExport={};
guidata(hObject, handles);



function Pindex_Callback(hObject, eventdata, handles)
% hObject    handle to Pindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pindex as text
%        str2double(get(hObject,'String')) returns contents of Pindex as a double


% --- Executes during object creation, after setting all properties.
function Pindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadProfile.
function LoadProfile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
index=str2num(get(handles.Pindex,'String'));
handles.Xall{index}=handles.x;
handles.Iall{index,1}=handles.smear_filt;
handles.Iall{index,2}=handles.smear;
set(handles.Pindex,'String',num2str(index+1));

AreaCheck=get(handles.AreaCheck,'Value');

if AreaCheck==1
    
    handles.AreasAll(index,:)=handles.Areas;
    handles.AreasNormAll(index,:)=handles.AreasNorm;
    set(handles.AreaTable,'Data',handles.AreasAll)
end

guidata(hObject, handles);


% --- Executes on button press in BGTemplate.
function BGTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to BGTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BGTemplate






% --- Executes during object creation, after setting all properties.
function BGTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BGTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.BGTemp=0;


% --- Executes on button press in mfilter.
function mfilter_Callback(hObject, eventdata, handles)
% hObject    handle to mfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
J=str2num(get(handles.mfparam,'string'));
handles.Gel_filt=medfilt2(handles.Gel,[J,J]);
axes(handles.image);
imshow(handles.Gel_filt);

guidata(hObject, handles);

function mfparam_Callback(hObject, eventdata, handles)
% hObject    handle to mfparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mfparam as text
%        str2double(get(hObject,'String')) returns contents of mfparam as a double


% --- Executes during object creation, after setting all properties.
function mfparam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mfparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Thresh_Median_Callback(hObject, eventdata, handles)
% hObject    handle to Thresh_Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Thresh_Median as text
%        str2double(get(hObject,'String')) returns contents of Thresh_Median as a double


% --- Executes during object creation, after setting all properties.
function Thresh_Median_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thresh_Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LeftBGButton.
function LeftBGButton_Callback(hObject, eventdata, handles)
% hObject    handle to LeftBGButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
It=handles.smear_filt;
S=str2num(get(handles.SmoothBox,'String'));
Start=str2num(get(handles.SmoothStart,'String'));
It(Start:end)=smooth(It(Start:end),S);
handles.LeftBG=It;
axes(handles.Rplot)
plot(handles.x,handles.LeftBG);
guidata(hObject, handles);
% --- Executes on button press in RightBGButton.
function RightBGButton_Callback(hObject, eventdata, handles)
% hObject    handle to RightBGButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
It=handles.smear_filt;
S=str2num(get(handles.SmoothBox,'String'));
Start=str2num(get(handles.SmoothStart,'String'));
It(Start:end)=smooth(It(Start:end),S);
handles.RightBG=It;
axes(handles.Rplot)
plot(handles.x,handles.RightBG);
guidata(hObject, handles);

% --- Executes on button press in SubtractBGProfile.
function SubtractBGProfile_Callback(hObject, eventdata, handles)
% hObject    handle to SubtractBGProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
Profile=(handles.RightBG+handles.LeftBG)/2;
Profile=Profile./Profile(1)*handles.smear_filt(1);
handles.smear_filt=handles.smear_filt-Profile;
axes(handles.Rplot);
plot(handles.smear_filt);
guidata(hObject, handles);

function SmoothBox_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmoothBox as text
%        str2double(get(hObject,'String')) returns contents of SmoothBox as a double


% --- Executes during object creation, after setting all properties.
function SmoothBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmoothBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SmoothStart_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SmoothStart as text
%        str2double(get(hObject,'String')) returns contents of SmoothStart as a double


% --- Executes during object creation, after setting all properties.
function SmoothStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SmoothStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in MedianFilter.
function MedianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MedianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MedianFilter


% --- Executes on button press in GetArea.
function GetArea_Callback(hObject, eventdata, handles)
% hObject    handle to GetArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AreaCheck.
function AreaCheck_Callback(hObject, eventdata, handles)
% hObject    handle to AreaCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AreaCheck


% --- Executes on button press in SumIntensity.
function SumIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to SumIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SumIntensity
check=get(handles.SumIntensity,'Value');
if check==1
    set(handles.MedianFilter,'Value',0);
end
guidata(hObject, handles);


% --- Executes on button press in ClearIntensity.
function ClearIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to ClearIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.Xall={};
handles.Iall={};
handles.AreasAll=[];
handles.AreasNormAll=[];
set(handles.Pindex,'String','1');
guidata(hObject, handles);
