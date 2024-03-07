% Interactively and visually set a threshold on a gray scale image.
% https://www.mathworks.com/matlabcentral/fileexchange/29372-thresholding-an-image?s_tid=srchtitle
%  File and function name : threshold
%  Written by    :  Image Analyst
%  Date : August 2006 - March 2022
%  Description   :  This program takes a color or monochrome image and lets
%                   the user threshold the monochrome image, or a single
%                   color band of a color image, via sliders to set the
%                   maximum and minimum thresholds.  The pixels in the
%                   thresholded range are shown in the middle image as a
%                   binary image (black/white), and the original image
%                   pixels are shown masked in the left image (gray scale
%                   or color).  The binary or masked image can be saved.
%                   Inputs are the low and high thresholds to start with
%                   and the image filename.
%                   The program returns the threshold values and the last
%                   color band that was used to select the threshold.
%
%  Program files : threshold.m, threshold.fig
%
%  Sample calls:
%	grayScaleImage = imread('coins.png');
%	[lowThreshold, highThreshold] = threshold(83, 255, grayScaleImage);
%	binaryImage = grayImage >= lowThreshold & grayImage <= highThreshold;
%
% or
%>>[lowThreshold, highThreshold, lastThresholdedBand] = threshold(10, 128,'D:\Images\MISC\MacbethSG.JPG'); % Read in color image.
%
% or
%>> monoImageArray = imread('coins.png');
%>>[lowThreshold, highThreshold, lastThresholdedBand] = threshold(83, 255, monoImageArray);
%
% or
%>> binsToSuppress = [71, 72]
%>>[lowThreshold, highThreshold, lastThresholdedBand] = threshold(83, 255, monoImageArray, binsToSuppress);
%
% or
% This version has the value, threshold2Use, that will determine if the min
% and max thresholds will be active (allowed to change).  If both are set to 0 there is an error.
%>> binsToSuppress = [71, 72]
%>> threshold2Use = [1, 0]
%>>[lowThreshold, highThreshold, lastThresholdedBand] = threshold(83, 255, monoImageArray, binsToSuppress, threshold2Use);
%
% or
% This version supplies an RGB image to use in display on the left hand side, instead of the gray level image that is normally there.
%>> binsToSuppress = [71, 72]
%>> threshold2Use = [1, 1]
%>>[lowThreshold, highThreshold, lastThresholdedBand] = threshold(83, 255, monoImageArray, binsToSuppress, threshold2Use, rgbImage);
%
% Example calling for a double, floating point monochrome image.
% The image does not need to be in the range 0-1 - it can have any range.
% Starting range is initialized to (-0.5 to 0.27).
%>> [lowThreshold highThreshold] = threshold(-0.5, 0.27, doubleImage);

function varargout = threshold(varargin)
% threshold Application M-file for threshold.fig
%    FIG = FIGMAINWINDOW launch threshold GUI.
%    FIGMAINWINDOW('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 03-Mar-2022 12:22:30
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @threshold_OpeningFcn, ...
	'gui_OutputFcn',  @threshold_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function threshold_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to threshold (see VARARGIN)
%clc;    % Clear Command Window.
global Original_Image;
global image_filename;
global axes1_image; % Can be a monochrome or color image array of any data type.
global axes2_image; % A uint8 image array.
global binsToSuppress;
global g_MaxThreshold;
global g_MinThreshold;
global numberOfColorBands;
global g_Floating;
global g_useMinThreshold;
global g_useMaxThreshold;
%     disp('Running threshold_OpeningFcn');
g_useMinThreshold = 1;
g_useMaxThreshold = 1;

%     g_MinThreshold = varargin(1);
%     g_MaxThreshold = varargin(2);
%     ImageFilename = varargin{3};
try
	% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
	warning('off', 'Images:initSize:adjustingMag');
	imageToDisplay = [];	% Default to not displaying a separately supplied image.
	
	numberOfArguments = nargin;
	if numberOfArguments <= 3
		% They ran this m-file with no arguments being passed in.
		errorMessage = sprintf('You must call the threshold() function by passing it an image array, a filename, and optionally some threshold values.');
		uiwait(warndlg(errorMessage));
		delete(handles.figMainWindow);
		return;
	elseif numberOfArguments == 4
		a = varargin{1};
		ImageFilename = a;
		if isa(a, 'char')
			% They passed in just the filename, but no starting thresholds.
			g_MinThreshold = [];
			g_MaxThreshold = [];
		elseif isnumeric(a)
			% They passed in the actual array, but no starting thresholds.
			g_MaxThreshold = max(max(a));
			g_MinThreshold = g_MaxThreshold/2;
		end
		
	elseif numberOfArguments <= 6
		% 3 input arguments will give nargin = 6 for some reason.
		[g_MinThreshold, g_MaxThreshold, ImageFilename] = deal(varargin{:}); % copies the contents of the cell array X to the separate variables Y1, Y2, Y3, ...
		binsToSuppress = [];
	elseif numberOfArguments == 7
		% Has an extra argument for the bins to supress.
		[g_MinThreshold, g_MaxThreshold, ImageFilename, binsToSuppress] = deal(varargin{:}); % copies the contents of the cell array X to the separate variables Y1, Y2, Y3, ...
	elseif numberOfArguments == 8
		% Has an extra argument for the bins to supress and an array of 2
		% values that determine if the min and max threshold is used,
		% respectfully
		[g_MinThreshold, g_MaxThreshold, ImageFilename, binsToSuppress, threshold2Use] = deal(varargin{:}); % copies the contents of the cell array X to the separate variables Y1, Y2, Y3, ...
		assert( numel(threshold2Use) == 2, 'The threshold2Use must be a 1x2 array.');
		assert( sum(threshold2Use) > 0, 'One of the threshold2Use needs to be set.');
		g_useMinThreshold = threshold2Use(1);
		g_useMaxThreshold = threshold2Use(2);
	elseif numberOfArguments == 9
		% They've passed in an image to display.  For example they want to display a color image imageToDisplay
		% but want to threshold on its gray scale version in ImageFilename.
		[g_MinThreshold, g_MaxThreshold, ImageFilename, binsToSuppress, threshold2Use, imageToDisplay] = deal(varargin{:}); % copies the contents of the cell array X to the separate variables Y1, Y2, Y3, ...
	end
	
	Original_Image =[];
	image_filename ='';
	axes1_image = [];
	axes2_image = [];
	axes(handles.axes1);
	fontSize = 12;
	
	% Store input array or file into image array variable.
	if isa(ImageFilename, 'char')
		%disp 'ImageFilename is a char';
		if exist(ImageFilename, 'file')
			LoadFile(handles, ImageFilename);
		else
			errorMessage = sprintf('The image file does not exist on disk.\n\nYou passed in this text string for the full filename of the disk file:\n%s', ImageFilename);
			msgboxw(errorMessage);
			return;
		end
	else
		%disp 'ImageFilename is not a string';
		% ImageFilename is really a numerical array, not a filename string.
		imageDimensions = size(ImageFilename);
		if length(imageDimensions) == 3
			numberOfColorBands = imageDimensions(3);
			xlabel('Color Band', 'FontSize', fontSize);
			handles.txtColorBand.Visible = 'on';
			handles.listbox_DisplayColorBand.Visible = 'on';
		else
			numberOfColorBands = 1;
			xlabel('Monochrome image', 'FontSize', fontSize);
			handles.txtColorBand.Visible = 'off';
			handles.listbox_DisplayColorBand.Visible = 'off';
			handles.axesHist.Visible = 'on';	% Make it visible (it was off for startup until and image is displayed.)
		end
		
		Original_Image = ImageFilename;
		axes1_image = Original_Image;
		axes(handles.axes1)
		warning off;
		if ~isempty(imageToDisplay)
			% Display the alternate image.
			imshow(imageToDisplay, []);
		else
			% Display the image indicated by ImageFilename.
			imshow(Original_Image, []);
		end
		SetUpImpixelinfo(handles);
		warning on;
	end
	% imshow() seems to take away the white color of the axes.  Restore it to white.
	ax = gca;
	ax.XColor = 'w';
	ax.YColor = 'w';

	% Determine if image is floating or integer
	if isa(Original_Image, 'double') || isa(Original_Image, 'single') || isa(Original_Image, 'float')
		g_Floating = 1;
	else
		g_Floating = 0;
	end
	
	% turn off the min threshold if need be
	if g_useMinThreshold == 0
		handles.slider_Min.Visible = 'off';
		handles.edit_Min.Visible = 'off';
		handles.text_Min.Visible = 'off';
	end
	
	% turn off the min threshold if need be
	if g_useMaxThreshold == 0
		handles.slider_Max.Visible = 'off';
		handles.edit_Max.Visible = 'off';
		handles.text_Max.Visible = 'off';
	end
	
	% Set up the slider range.
	maxSliderValue = max(max(Original_Image));
	minSliderValue = min(min(Original_Image));
	if maxSliderValue == minSliderValue
		% If the image has no range, it's just a uniform gray level.
		maxSliderValue = minSliderValue + 1;
	end
	
	% The above two will be 1x1x3 arrays if the image is a color image.
	% Set sliders to the extreme if it's color.
	if numberOfColorBands > 1
		maxSliderValue = max(maxSliderValue);
		minSliderValue = min(minSliderValue);
	end
	
	range = maxSliderValue - minSliderValue;
	% Specify the amount to move when they click the scroll bar arrow.
	% Determine if image is floating or integer
	if g_Floating
		slider_step(1) = 0.01;
		slider_step(2) = 0.1;
	else
		if isa(Original_Image, 'uint8')
			% Set up 8 bit images to have arrow move by 1 gray level, and clicking in slider moving it by 10 gray levels.
			slider_step(1) = 1.0 / double(range);
			slider_step(2) = 10.0 * slider_step(1);
		else
			% Set up 16 bit images to have arrow move by 10 gray level, and clicking in slider moving it by 1000 gray levels.
			slider_step(1) = 10.0 / double(range);
			slider_step(2) = 100.0 * slider_step(1);
		end
	end
	if slider_step(2) > 1
		% Can happen is the range is small
		slider_step(1) = 0.01;
		slider_step(2) = 0.1;
	end
	
	% Warning: slider control can not have a Value outside of Min/Max range
	% Control will not be rendered until all of its parameter values are valid.
	% also, if either the min or max is not being used, set it to the
	% limits of the image.
	if g_MinThreshold < minSliderValue || g_MinThreshold > maxSliderValue || g_useMinThreshold == 0
		g_MinThreshold = minSliderValue;
	end
	
	if g_MaxThreshold < minSliderValue || g_MaxThreshold > maxSliderValue || g_useMaxThreshold == 0
		g_MaxThreshold = maxSliderValue;
	end
	
	% If they're the same, make the max of the range be the max available.
	if g_MaxThreshold == g_MinThreshold
		g_MaxThreshold = maxSliderValue;	% Select bright objects.
	end
	
	handles.slider_Min.SliderStep = slider_step;
	handles.slider_Min.Max = maxSliderValue;
	handles.slider_Min.Min = minSliderValue;
	handles.slider_Min.Value = g_MinThreshold;
	handles.slider_Max.SliderStep = slider_step;
	handles.slider_Max.Max = maxSliderValue;
	handles.slider_Max.Min = minSliderValue;
	handles.slider_Max.Value = g_MaxThreshold;
	handles.edit_Min.String = g_MinThreshold;
	handles.edit_Max.String = g_MaxThreshold;
	
	% Plot the histogram
	PlotHistogram(handles, Original_Image, binsToSuppress);
	
	ShowThresholdedBinaryImage(hObject, eventdata, handles);
	% Enlarge figure to full screen.
	g = gcf;
	g.WindowState = 'maximized';
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
% Call guidata anytime you use the set() command or directly access
% a variable you added to the handles structure.  Generally it's a good
% practice to just automatically add this at the end of every function.
guidata(hObject, handles);

% Wait for user to interact with this macro.  Macro will continue
% with the statement after uiwait when the user exits the GUI.
uiwait();
%%%% Continues here after user exits GUI (and a uiresume is issued)
% For example, if the btnClose button's callback has these lines.
% uiresume(handles.figMainWindow);
% delete(handles.figMainWindow);
% Then it will resume from here and run run threshold_OutputFcn before shutting down.
% If some other function calls uiresume() without deleting the handle, then
% it will resume here, and call the threshold_OutputFcn but will stay up and running.

% End of threshold_OpeningFcn function
%=====================================================================================================================================================

% --- Outputs from this function are returned to the command line.
function varargout = threshold_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global g_MaxThreshold;
global g_MinThreshold;
global g_LastThresholdedColorBand;
varargout{1} = g_MinThreshold;
varargout{2} = g_MaxThreshold;
varargout{3} = g_LastThresholdedColorBand;
% handles.figMainWindow.Color = [0, 0, .25];
% g = gcf;
% g.WindowState = 'maximized';
%disp('Exited threshold macro');

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and
%| sets objects' callback properties to call them through the FEVAL
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'threshold_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figMainWindow, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


%=====================================================================================================================================================
function slider_Max_Callback(hObject, eventdata, handles)
% Stub for Callback of the uicontrol handles.slider_Max.
global g_Floating;
global g_MaxThreshold;
global g_MinThreshold;
try
	% Get the value of the sliders
	minSliderValue = handles.slider_Min.Value;
	maxSliderValue = handles.slider_Max.Value;
	% Round it to an integer if the image is a uint (not a floating).
	if g_Floating == 0
		maxSliderValue = round(maxSliderValue);
	end
	
	%disp(['Moved max slider.  New value = ' num2str(maxSliderValue)]);
	
	% Make sure it's not greater than the Max slider.
	if (maxSliderValue < g_MinThreshold)
		% Illegal value. Clip to maxSliderValue
		maxSliderValue = g_MinThreshold;
		%disp(['Adjusted max slider.  New value = ' num2str(maxSliderValue)]);
	end
	% Set the global variable.
	g_MaxThreshold = maxSliderValue; % Update global variable.
	handles.slider_Max.Value = maxSliderValue;
	
	% Set up both edit fields according to what their slider's value is.
	if rem(minSliderValue, 1) == 0
		% It's an integer
		handles.edit_Min.String = sprintf('%d', minSliderValue);
	else
		% It's floating point.
		handles.edit_Min.String = sprintf('%.4f', minSliderValue);
	end
	if rem(maxSliderValue, 1) == 0
		% It's an integer
		handles.edit_Max.String = sprintf('%d', maxSliderValue);
	else
		% It's floating point.
		handles.edit_Max.String = sprintf('%.4f', maxSliderValue);
	end

	% Call guidata anytime you use the set() command or directly access a variable you added to the handles structure.  
	% Generally it's a good practice to just automatically add this at the end of every function.
	guidata(hObject, handles);
	
	% Apply the threshold to the display so we can see its effect.
	ShowThresholdedBinaryImage(hObject, eventdata, handles);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from slider_Max_Callback()


%=====================================================================================================================================================
% --- Executes on slider movement.
function slider_Min_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global g_Floating;
global g_MaxThreshold;
global g_MinThreshold;
try
	% Get the value of the sliders
	minSliderValue = handles.slider_Min.Value;
	maxSliderValue = handles.slider_Max.Value;
	% Round it to an integer if the image is a uint (not a floating).
	if g_Floating == 0
		minSliderValue = round(minSliderValue);
	end
	
	%disp(['Moved min slider.  New value = ' num2str(minSliderValue)]);
	
	% Make sure it's not greater than the Max slider.
	if (minSliderValue > g_MaxThreshold)
		% Illegal value. Clip to maxSliderValue
		minSliderValue = g_MaxThreshold;
		%disp(['Adjusted min slider.  New value = ' num2str(minSliderValue)]);
	end
	handles.slider_Min.Value = minSliderValue;
	% Set the global variable.
	g_MinThreshold = minSliderValue;

	% Set up both edit fields according to what their slider's value is.
	if rem(minSliderValue, 1) == 0
		% It's an integer
		handles.edit_Min.String = sprintf('%d', minSliderValue);
	else
		% It's floating point.
		handles.edit_Min.String = sprintf('%.4f', minSliderValue);
	end
	if rem(maxSliderValue, 1) == 0
		% It's an integer
		handles.edit_Max.String = sprintf('%d', maxSliderValue);
	else
		% It's floating point.
		handles.edit_Max.String = sprintf('%.4f', maxSliderValue);
	end

	% Call guidata anytime you use the set() command or directly access
	% a variable you added to the handles structure.  Generally it's a good
	% practice to just automatically add this at the end of every function.
	guidata(hObject, handles);
	
	% Apply the threshold to the display so we can see its effect.
	ShowThresholdedBinaryImage(hObject, eventdata, handles);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from slider_Min_Callback()


%=====================================================================================================================================================
function edit_Max_Callback(hObject, eventdata, handles)
% Stub for Callback of the uicontrol handles.edit_Max.
global g_MaxThreshold;
global g_MinThreshold;
global g_LastThresholdedColorBand;
try
	g_LastThresholdedColorBand = handles.listbox_DisplayColorBand.Value - 1;
	
	% Get the value from the control.
	max_value = str2double(handles.edit_Max.String);
	if (max_value <= 255) && (max_value >= g_MinThreshold) && isnan(max_value) == false
% 		g_MaxThreshold = round(max_value); % Only for integer data types, if that.
		g_MaxThreshold = max_value;
		handles.slider_Max.Value = g_MaxThreshold;
	end
	handles.edit_Max.String = g_MaxThreshold;
	% Call guidata anytime you use the set() command or directly access
	% a variable you added to the handles structure.  Generally it's a good
	% practice to just automatically add this at the end of every function.
	guidata(hObject, handles);
	ShowThresholdedBinaryImage(hObject, eventdata, handles);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from edit_Max_Callback()


%=====================================================================================================================================================
function edit_Min_Callback(hObject, eventdata, handles)
% Stub for Callback of the uicontrol handles.edit_Min.
global g_MaxThreshold;
global g_MinThreshold;
global g_LastThresholdedColorBand;
try
	g_LastThresholdedColorBand = handles.listbox_DisplayColorBand.Value - 1;
	
	min_value = str2double(handles.edit_Min.String);
	if (min_value >= 0) && (min_value <= g_MaxThreshold) && isnan(min_value) == false
% 		g_MinThreshold = round(min_value); % Only for integer data types, if that.
		g_MinThreshold = min_value;
		handles.slider_Min.Value = g_MinThreshold;
	end
	handles.edit_Min.String = g_MinThreshold;
	% Call guidata anytime you use the set() command or directly access
	% a variable you added to the handles structure.  Generally it's a good
	% practice to just automatically add this at the end of every function.
	guidata(hObject, handles);
	ShowThresholdedBinaryImage(hObject, eventdata, handles);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from edit_Min_Callback()


%=====================================================================================================================================================
% Read in a file from disk.
function LoadFile(handles, filename)
% Stub for Callback of the uicontrol handles.pushbutton_LoadFile.
global Original_Image;
global image_filename;
global axes1_image; % Can be a monochrome or color image array of any data type.
global numberOfColorBands;
try
	if (strcmp(filename,image_filename) == 0)
		%reply = questdlg('Sure that u want to load this file?','Loading File...','Yes','No','No');
		reply = 'Yes';  % Don't bother asking - just go ahead and do it.
		if (strcmp(reply,'Yes') == 1)
			try %Checks the properties of the file
				info = imfinfo(filename);
				% Disable scroll bars.
				handles.slider_Min.Enable = 'off';
				handles.slider_Max.Enable = 'off';
				handles.edit_Min.Enable = 'off';
				handles.edit_Max.Enable = 'off';
				numberOfColorBands = info.BitDepth / 8;
				switch info.ColorType
					case 'truecolor'
						errormessage = '';
						handles.listbox_DisplayColorBand.Visible = 'on';
						handles.txtColorBand.Visible = 'on';
						% Turn off the other two images until they click a color
						% band.
						handles.axes2.Visible = 'off';
						handles.axes3.Visible = 'off';
					case 'grayscale'
						errormessage = '';
						handles.listbox_DisplayColorBand.Visible = 'off';
						xlabel('Monochrome image', 'FontSize', 12);
						handles.txtColorBand.Visible = 'off';
						% Enable scroll bars.
						handles.slider_Min.Enable = 'on';
						handles.slider_Max.Enable = 'on';
						handles.edit_Min.Enable = 'on';
						handles.edit_Max.Enable = 'on';
					case 'indexed'
						errormessage = 'MATLAB reports that this is an INDEXED image. ThProgram accepts TRUECOLOR only';
					otherwise
						errormessage = 'Unknown image color system';
				end
			catch ME
				errormessage = 'This is NOT an IMAGE file';
			end
			if isempty(errormessage)
				try
					Original_Image = imread(filename);
					axes(handles.axes1)
					warning off;
					imshow(Original_Image);
					SetUpImpixelinfo(handles);
					% imshow() seems to take away the white color of the axes.  Restore it to white.
					ax = gca;
					ax.XColor = 'w';
					ax.YColor = 'w';
					warning on;
					axes1_image = Original_Image;
					image_filename = filename;
				catch ME
					errorMessage = sprintf('Error reading and displaying image file\.Error in function %s() at line %d.\n\nError Message:\n%s', ...
						ME.stack(1).name, ME.stack(1).line, ME.message);
					WarnUser(errorMessage);
				end
			else
				msgboxw(errormessage);
			end
		end
	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from LoadFile()


%=====================================================================================================================================================
% --- Executes on selection change in listbox_DisplayColorBand.
function listbox_DisplayColorBand_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_DisplayColorBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_DisplayColorBand contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_DisplayColorBand
global Original_Image;
global axes1_image; % Can be a monochrome or color image array of any data type.
try
	fontSize = 12;
	% Now that they've clicked on something, show the other two images.
	handles.axes2.Visible = 'on';
	handles.axes3.Visible = 'on';
	handles.txtInfo.String = '';   % Clear the info text.
	% Find out what color band they've requested.
	colorBand = handles.listbox_DisplayColorBand.Value;
	axes(handles.axes1);    % Switch current figure to axes2, the middle picture.
	errormessage = '';
	if colorBand == 1
		% Disable scroll bars.
		handles.slider_Min.Enable = 'off';
		handles.slider_Max.Enable = 'off';
		handles.edit_Min.Enable = 'off';
		handles.edit_Max.Enable = 'off';
	else
		% Enable scroll bars.
		handles.slider_Min.Enable = 'on';
		handles.slider_Max.Enable = 'on';
		handles.edit_Min.Enable = 'on';
		handles.edit_Max.Enable = 'on';
		handles.axesHist.Visible = 'on';	% Make it visible (it was off for startup until and image is displayed.)
	end
	
	switch colorBand
		case 1 % Full color
			axes1_image = Original_Image;
			xlabel('Full color', 'FontSize', fontSize);
		case 2 % Red Band
			axes1_image = Original_Image(:,:,1);
			xlabel('Red Band', 'FontSize', fontSize);
		case 3 % Green Band
			axes1_image = Original_Image(:,:,2);
			xlabel('Green Band', 'FontSize', fontSize);
		case 4 % Blue Band
			axes1_image = Original_Image(:,:,3);
			xlabel('Blue Band', 'FontSize', fontSize);
		otherwise
			errormessage = 'Invalid selection';
	end
	if colorBand >= 2 && colorBand <= 4
		% Plot the histogram
		PlotHistogram(handles, axes1_image, []);
		ShowThresholdedBinaryImage(hObject, eventdata, handles);
	end
	if isempty(errormessage)
		warning off;
		imshow(axes1_image);
		% imshow() seems to take away the white color of the axes.  Restore it to white.
		ax = gca;
		ax.XColor = 'w';
		ax.YColor = 'w';
		SetUpImpixelinfo(handles);
		warning on;
	else
		xlabel('Error displaying');
	end
catch ME
	errorMessage = sprintf('Error viewing selected color band\nError in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	handles.txtInfo.String = errormessage;
	WarnUser(errorMessage);
end
% Call guidata anytime you use the set() command or directly access
% a variable you added to the handles structure.  Generally it's a good
% practice to just automatically add this at the end of every function.
guidata(hObject, handles);
ShowThresholdedBinaryImage(hObject, eventdata, handles);


%=====================================================================================================================================================
% Performs the threshold and shows the binary, thresholded image plus the masked image.
function ShowThresholdedBinaryImage(hObject, eventdata, handles)
global g_MaxThreshold;
global g_MinThreshold;
global axes1_image; % Can be a monochrome or color image array of any data type.
global axes2_image; % A uint8 image array.
global g_LastThresholdedColorBand;
global numberOfColorBands;
global axesShowingInMiddleImage;
try
	handles = guidata(hObject);
	fontSize = 12;
	axes(handles.axes2);
	colorBand = handles.listbox_DisplayColorBand.Value;
	if (colorBand ~= 1) || (numberOfColorBands == 1)
		% If they're displaying one color band, or a monochrome image...
		handles.txtInfo.String = 'Finding thresholding regions...';
		% Set all pixels which are in the threshold selection range to white/true/1,
		% and all pixels which are outside the threshold selection range to black/false/0.
		axes2_image = (g_MinThreshold <= axes1_image) & (axes1_image <= g_MaxThreshold);
		% axes2_image will be a binary (logical) image with values of 0 and/or 1 (false or true) only.
		
		message = ['Masked image between ',num2str(g_MinThreshold),' and  ',...
			num2str(g_MaxThreshold)];
		xlabel(handles.axes2, message, "FontSize", fontSize);
		try
			handles.txtInfo.String = 'Displaying Axes2 image';
			imshow(axes2_image, [], 'Parent', handles.axes2);
			% imshow() seems to take away the white color of the axes.  Restore it to white.
			ax = gca;
			ax.XColor = 'w';
			ax.YColor = 'w';
			handles.txtInfo.String = '';
			SetUpImpixelinfo(handles);
		catch ME
			handles.txtInfo.String = 'Error displaying Axes2';
			errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
				ME.stack(1).name, ME.stack(1).line, ME.message);
			WarnUser(errorMessage);
		end
		g_LastThresholdedColorBand = handles.listbox_DisplayColorBand.Value - 1;
		PlaceThresholdBar(handles, g_MinThreshold, g_MaxThreshold);
	end
	axesShowingInMiddleImage = 2;	% Flag that we're showing the binary image.
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
% Call guidata anytime you use the set() command or directly access
% a variable you added to the handles structure.  Generally it's a good
% practice to just automatically add this at the end of every function.
guidata(hObject, handles);

ShowMaskedImage(hObject, eventdata, handles);
guidata(hObject, handles);
return; % ShowThresholdedBinaryImage


%=====================================================================================================================================================
function handles = ShowMaskedImage(hObject, eventdata, handles)
global axes1_image; % Can be a monochrome or color image array of any data type.
global axes2_image; % A uint8 image array.
global maskedImage;
global numberOfColorBands;
try
	radioButtonSetting = handles.radBlackOutside.Value;
	
	% You can only multiply integers if they are of the same type.
	% Convert the type of axes2_image (the binary mask) to the same data type as axes1_image (the original image).
	axes2_image = cast(axes2_image, class(axes1_image));  % Convert type of axes2_image to same type as axes1_image.
	
	handles.txtInfo.String = '';
	maxSliderValue = handles.slider_Max.Max;
	minSliderValue = handles.slider_Min.Min;
	if (~isempty(axes1_image) && ~isempty(axes2_image))
		handles.txtInfo.String = 'Displaying images';
		colorBand = handles.listbox_DisplayColorBand.Value;
		if radioButtonSetting == 1
			% Make image black outside the threshold range.
			if (colorBand ~= 1) || (numberOfColorBands == 1)
				% Single color band.
				maskedImage = axes2_image .* axes1_image;    % maskedImage has the same data type as axes1_image.
			else
				% Full color image.
				maskedImage(:,:,1) = uint8(axes1_image(:,:,1) .* axes2_image);
				maskedImage(:,:,2) = uint8(axes1_image(:,:,2) .* axes2_image);
				maskedImage(:,:,3) = uint8(axes1_image(:,:,3) .* axes2_image);
			end
		else
			% Make image white outside the threshold range.
			maxValueForThisDataType = maxSliderValue;
			% 			maxValueForThisDataType = intmax(strDataType);
			dblInvertedBinary = maxValueForThisDataType * (1 - axes2_image); % Note: axes2_image is still in the range 0-1.
			% You can only multiply integers if they are of the same type.
			% Convert the type of axes2_image (the binary mask) to the same data type as axes1_image (the original image).
			dblInvertedBinary = cast(dblInvertedBinary, class(axes1_image));  % Convert type of dblInvertedBinary to same type as axes1_image.
			
			if (colorBand ~= 1) || (numberOfColorBands == 1)
				% Single color band.
				maskedImage = axes2_image .* axes1_image + dblInvertedBinary;
			else
				% Full color image.
				maskedImage(:,:,1) = uint8(axes1_image(:,:,1) .* axes2_image + dblInvertedBinary);
				maskedImage(:,:,2) = uint8(axes1_image(:,:,2) .* axes2_image + dblInvertedBinary);
				maskedImage(:,:,3) = uint8(axes1_image(:,:,3) .* axes2_image + dblInvertedBinary);
			end
		end
		
		try
			% Display the image.
			% First see if they want to scale the intensity to make it brighter.
			chkScaleIntensity = handles.chkScaleIntensity.Value;
			warning off;
			axes(handles.axes3);
			if chkScaleIntensity && (numberOfColorBands == 1)
				imshow(maskedImage, []);
			else
				if isa(maskedImage, 'single') || isa(maskedImage, 'double')
					scaledImage = im2single(maskedImage);
					imshow(maskedImage, [minSliderValue maxSliderValue]);
				else
					imshow(maskedImage);
				end
			end
			SetUpImpixelinfo(handles);
			% imshow() seems to take away the white color of the axes.  Restore it to white.
			ax = gca;
			ax.XColor = 'w';
			ax.YColor = 'w';
			warning on;
			handles.txtInfo.String = '';
		catch ME
			handles.txtInfo.String = 'Error displaying masked image.';
			errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
				ME.stack(1).name, ME.stack(1).line, ME.message);
			WarnUser(errorMessage);
		end
	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
guidata(hObject, handles);
return; % from ShowMaskedImage()


%=====================================================================================================================================================
function PlotHistogram(handles, imageArray, binsToSuppress)
global g_Floating;
try
	imageDimensions = size(imageArray);
	if length(imageDimensions) == 3
		% It's a color imageArray.  We can't plot it.
		handles.axesHist.Visible = 'off';	% Make it invisible.
	else
		% It's a monochrome imageArray.  We can plot it.
		minSliderValue = handles.slider_Min, 'min';
		maxSliderValue = handles.slider_Min, 'max';
		% Find out how many bins we need for this data type so that we get one gray level per bin.
		strDataType = class(imageArray); % Get data type of original_image.
		minGrayLevelInImage = min(min(imageArray));
		maxGrayLevelInImage = max(max(imageArray));
		if g_Floating
			% It's a floating point image, not an integer image that is our "starting image."
			% However can't use imhist on abritrary floating point images.
			% They need to be in the range 0-1, so we need to scale to 0-1 to make sure it will work.
			rangeS = maxGrayLevelInImage - minGrayLevelInImage;
			normalizedImage = (imageArray - minGrayLevelInImage) / rangeS;
			% Check to verify that range is now 0-1.
			minValueNorm = min(min(normalizedImage));
			maxValueNorm = max(max(normalizedImage));
			
			% Let's get its histogram into 256 bins.
			[counts, grayLevelsD] = imhist(normalizedImage, 256);
			% But now grayLevelsD goes from 0 to 1.
			% We want it to go from the original range, so we need to scale.
			grayLevels = rangeS * grayLevelsD + minGrayLevelInImage;
		else
			% Uint8 or uint16
			if contains(strDataType, '8')
				numberOfBins = 256;
			else
				numberOfBins = 512;
			end
% 			maxValueForThisDataType = intmax(strDataType);
% 			minValueForThisDataType = intmin(strDataType);
% 			numberOfBins = double(maxValueForThisDataType) - double(minValueForThisDataType) + 1; % imhist() requires that this be a double data type.
			% Get a histogram of the entire image.
			[counts, grayLevels] = imhist(imageArray, numberOfBins);    % make sure you label the axes after imhist because imhist will destroy them.
		end
		
		% Find the last non-zero bin so we can plot just up to there to get better horizontal resolution.
		lastBinUsed = find(counts > 0, 1, 'last');
		% Find the first non-zero bin.
		firstBinUsed = find(counts >= 0, 1, 'first');
		% Get subportion of array that has non-zero data.
		counts = counts(firstBinUsed : lastBinUsed);
		grayLevels = grayLevels(firstBinUsed : lastBinUsed);
		
		% Suppress any bins in the histogram, if the user asked us to.
		if length(binsToSuppress)
			if g_Floating == 0
				% Integer gray levels.
				for k = 1: length(binsToSuppress)
					grayLevelToSuppress = binsToSuppress(k);
					% Find out what bin that gray level occurs at.
					bin = find(grayLevels == grayLevelToSuppress, 1, 'first');
					counts(bin) = 0;
				end
			else
				% Floating point image.  Just zero out the bins directly.
				counts(binsToSuppress) = 0;
			end
		end
		% Get rid of outliers.
		% Find the max of the two bins on either side of a bin
		maxOfNeighborBins = imdilate(counts, [1;0;1]);
		% Find bins that are more than 5 times the max of the neighbors.
		tallSpikes = counts > 5 * maxOfNeighborBins;
		% Set those bins to the max of the neighbor bins.
		counts(tallSpikes) = maxOfNeighborBins(tallSpikes);
		
		% Plot the histogram in the histogram viewport.
		handles.axesHist.Visible = 'on';	% Make it visible (it was off for startup until and image is displayed.)
		axes(handles.axesHist);  % makes existing axes handles.axesHist the current axes.
		% Plot the histogram as a line curve.
		if handles.radLogCount.Value
			semilogy(grayLevels, counts);
		else
			plot(grayLevels, counts, 'LineWidth', 2);
		end

		% Set up the histogram axes
		axHist = handles.axesHist;

		% Plot grid lines if requested.
		showGridLines = handles.chkGridLines.Value;
		if showGridLines
			grid(axHist, 'on');
		else
			grid(axHist, 'off');
		end

		% Change the colors of the histogram axes.
		axHist.XColor = 'g';
		xlabel('Gray Level', 'Color', 'white');
		axHist.YColor = 'g';
		ylabel('Pixel Count', 'Color', 'white');

		% Make the opacity of the grid lines more so they show up better.  Right now the alpha is only 0.15 so they are faint.
		axHist.GridAlpha = 0.3;
		axHist.GridColor = 'k';
		% Set up the limits
		xlim([0, maxGrayLevelInImage]);
% 		if ~g_Floating
% 			xlim([0 maxValueForThisDataType]);
% 		end

		% Set up tick marks on the histogram axes.  Originally they are so wide apart.
		if contains(class(imageArray), 'uint8')
			xticks(axHist, [0:20:255, 255]);
		elseif contains(class(imageArray), 'uint16')
			xticks(axHist, [0:2000:65535, 65535]);
		elseif contains(class(imageArray), 'double')
			xt = linspace(grayLevels(1), grayLevels(end), 100);
			xticks(axHist, xt);
		end

		% Set up the axes tick mark colors on the three images to be white.
		ax = handles.axes1;
		ax.XColor = 'w';
		ax.YColor = 'w';
		ax = handles.axes2;
		ax.XColor = 'w';
		ax.YColor = 'w';
		ax = handles.axes3;
		ax.XColor = 'w';
		ax.YColor = 'w';

	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return;  % End of PlotHistogram()


%=====================================================================================================================================================
% Deletes any existing threshold bar from the histogram plot.
% Plots the new bar at the location of dblNewThreshold
function PlaceThresholdBar(handles, dblLowThreshold, dblHighThreshold)
global g_useMinThreshold;
global g_useMaxThreshold;
global hNewLowBar;
global hNewHighBar;
try
	% If the histogram plot is visible, replot the threshold line
	if handles.axesHist.Visible
		axes(handles.axesHist);  % makes existing axes handles.axesHist the current axes.
		% Set up the horizontal axis to be 0 - 255 in range.
		XRange = xlim(handles.axesHist);
		%XRange(2) = 255;
		xlim(handles.axesHist, XRange);
		% Make sure threshold bar location is in the range
		% along the horizontal X (gray level) axis.
		XRange = xlim(handles.axesHist);
		barWidth = XRange(2) / 200;
		maxXValue = XRange(2);
		if dblLowThreshold > maxXValue
			dblLowThreshold = maxXValue;
		end
		if dblHighThreshold > maxXValue
			dblHighThreshold = maxXValue;
		end
		
		% Erase the old lines.
		% Need to specify handles in findobj, otherwise it will get hggroups in other figures,
		% and then the delete() function will blow them away, believe it or not!
		% 		hOldLowBar=findobj(handles.axesHist, 'type', 'hggroup');
		% 		delete(hOldLowBar);
		% 		hOldHighBar=findobj(handles.axesHist, 'type', 'hggroup');
		% 		delete(hOldHighBar);
		delete(hNewLowBar);
		delete(hNewHighBar);
		
		% Draw a vertical line at the threshold to show where we're splitting the histogram.
		hold on;
		YRange = ylim(handles.axesHist);
		maxYValue = YRange(2);
		% Draw a bar at the low threshold.
		if g_useMinThreshold
			hNewLowBar = xline(dblLowThreshold, 'Color', 'r', 'LineWidth', 2);
		end
		% Draw another bar at the high threshold.
		if g_useMaxThreshold
			hNewHighBar = xline(dblHighThreshold, 'Color', 'r', 'LineWidth', 2);
		end
		% WEIRD QUIRK: for some gray levels, the bar doesn't get drawn.
		% May be some kind of aliasing/subsampling problem.
		hold off;
	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return;	% End of PlaceThresholdBar()


%=====================================================================================================================================================
% Saves the masked (right hand) image.
function mnuFileSaveMasked_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFileSaveMasked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global maskedImage;
global image_filename;
try
	if ~isempty(maskedImage)
		% Get the folder where the original image resides.
		[folder, baseFileName, ext] = fileparts(image_filename);
		% Ask user to specify file name and folder.
		[fileName, folder] = uiputfile([folder '\MaskedImage.png'], 'Please specify the folder and name for the masked image');
		NameOfFileToSave = fullfile(folder, fileName);
		% Save the mask image for later recall if desired.
		imwrite(maskedImage, NameOfFileToSave);
		handles.txtInfo.String = ['Saved masked image ' NameOfFileToSave];
	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return;  % End of mnuFileSaveMasked_Callback()


%=====================================================================================================================================================
% Saves the binary (middle) image.
function mnuFileSaveBinary_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFileSaveBinary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global axes2_image; % A uint8 image array.
global image_filename;
try
	if ~isempty(axes2_image)
		% Get the folder where the original image resides.
		[folder, baseFileName, ext] = fileparts(image_filename);
		% Ask user to specify file name and folder.
		[fileName, folder] = uiputfile([folder '\BinaryImage.png'], 'Please specify the folder and name for the binary image');
		NameOfFileToSave = fullfile(folder, fileName);
		% Save the mask image for later recall if desired.
		imwrite(axes2_image, NameOfFileToSave);
		handles.txtInfo.String = ['Saved binary image ' NameOfFileToSave];
	end
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return;  % End of mnuFileSaveBinary_Callback()

%=====================================================================================================================================================
% --- Executes on button press in pbMeasure.
function pbMeasure_Callback(hObject, eventdata, handles)
% hObject    handle to pbMeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% You must declare global variables within functions
% in order for the functions to see them.
global axes1_image; % Can be a monochrome or color image array of any data type.
global axes2_image; % A uint8 image array.
try
	% Use bw to mask off the original.
	%maskedImage = uint8(double(axes2_image) .* double(axes1_image));
	maskedImage = axes1_image .* uint8(axes2_image);
	% Display masked image.
	%figure; % Create new popup window. Otherwise it show up on the main form.
	%imshow(maskedImage);
	% Get a column vector of all the pixels within the masked area.
	maskedPixels = maskedImage (logical(axes2_image));   % Gives a column vector.
	%uiwait(msgbox('Click OK to see masked original image.'));
	meanGL = mean(maskedPixels(:));
	strLine1 = ['The mean GL is ', num2str(meanGL)];
	numberOfThresholdedPixels = length(maskedPixels(:));
	strLine2 =['The number of thresholded pixels is ', num2str(numberOfThresholdedPixels)];
	[imageWidth, imageHeight] = size(axes1_image);
	numberOfOriginalPixels = imageWidth * imageHeight;
	strLine3 = ['The number of pixels in the original image is ', num2str(numberOfOriginalPixels)];
	areaFraction = double(numberOfThresholdedPixels) / double(numberOfOriginalPixels) ;
	strLine4 =['The area fraction is ', num2str(areaFraction)];
	caResults = {strLine2;strLine3;strLine4;strLine1}; % Create cell array to get on multiple lines.
	%msgboxw(caResults);
	handles.txtResults.String = caResults;
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return;  % End of pbMeasure_Callback()


%=====================================================================================================================================================
function msgboxw(in_strMessage)
uiwait(msgbox(in_strMessage));
return


%=====================================================================================================================================================
% --- Executes on button press in btnOK.
function btnOK_Callback(hObject, eventdata, handles)
% hObject    handle to btnOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figMainWindow);
delete(handles.figMainWindow);


%=====================================================================================================================================================
function mnuFile_Callback(hObject, eventdata, handles)


%=====================================================================================================================================================
function mnuFileExit_Callback(hObject, eventdata, handles)
uiresume(handles.figMainWindow);
delete(handles.figMainWindow);


%=====================================================================================================================================================
% --- Executes when user attempts to close threshold.
function figMainWindow_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure
delete(hObject);


%=====================================================================================================================================================
% --- Executes on button press in radWhiteOutside.
function radWhiteOutside_Callback(hObject, eventdata, handles)
% MATLAB quirk: when you repeatedly click the radio button, it toggles it's value rather than
% keeping it set (selected) as it should.  As a workaround, make sure it gets set (not toggled) and
% make sure the other radio buttons are set to false.
% Change mouse pointer (cursor) to an hourglass.
% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
g =gcf
g.Pointer = 'watch';
drawnow;	% Cursor won't change right away unless you do this.
try
	handles.radWhiteOutside.Value = true; % Set this one true.
	% Set the other radio button's value the opposite of this one.
	handles.radBlackOutside.Value = false;
	ShowMaskedImage(hObject, eventdata, handles);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
% Change mouse pointer (cursor) to an arrow.
g.Pointer = 'arrow';
drawnow;	% Cursor won't change right away unless you do this.
guidata(hObject, handles);
return;  % End of radWhiteOutside_Callback()


%=====================================================================================================================================================
% --- Executes on button press in radBlackOutside.
function radBlackOutside_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radBlackOutside
% MATLAB quirk: when you repeatedly click the radio button, it toggles it's value rather than
% keeping it set (selected) as it should.  As a workaround, make sure it gets set (not toggled) and
% make sure the other radio buttons are set to false.
% Change mouse pointer (cursor) to an hourglass.
% QUIRK: use 'watch' and you'll actually get an hourglass not a watch.
g = gcf;
g.Pointer = 'watch';
drawnow;	% Cursor won't change right away unless you do this.
try
	handles.radBlackOutside.Value = true; % Set this one true.
	% Set the other radio button's value the opposite of this one.
	handles.radWhiteOutside.Value = false;
	ShowMaskedImage(hObject, eventdata, handles);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
% Change mouse pointer (cursor) to an arrow.
g.Pointer = arrow';
drawnow;	% Cursor won't change right away unless you do this.
guidata(hObject, handles);
return;  % End of radBlackOutside_Callback()


%=====================================================================================================================================================
% --- Executes on button press in btnFlicker.
function btnFlicker_Callback(hObject, eventdata, handles)
global axes1_image; % Can be a monochrome or color image array of any data type.
global axesShowingInMiddleImage; % 1 = showing original image, 2  showing binary image.
try
	axes(handles.axes2);
	fontSize = 12;
	if axesShowingInMiddleImage == 1 % 1 = showing original image, 2  showing binary image.
		% The original image is currently showing.
		% Switch to showing the binary image.
		ShowThresholdedBinaryImage(hObject, eventdata, handles);
		axesShowingInMiddleImage = 2; % 1 = showing original image, 2  showing binary image.
		xlabel(handles.axes2, 'Binary Image', 'FontSize', fontSize, 'Color', 'w', 'FontWeight','bold');
	else
		% The binary image is currently showing.
		% Switch to showing the original image.
		imshow(axes1_image, []);
		SetUpImpixelinfo(handles);
		axesShowingInMiddleImage = 1; % 1 = showing original image, 2  showing binary image.
		% 		axes1String = handles.textColorBand.String;
		xlabel(handles.axes2, 'Original Image', 'FontSize', fontSize, 'Color', 'w', 'FontWeight','bold');
	end
	% imshow() seems to take away the white color of the axes.  Restore it to white.
	ax = handles.axes2;
	ax.XColor = 'w';
	ax.YColor = 'w';
	xlabel(handles.axes1, 'Original Image', 'FontSize', fontSize, 'Color', 'w', 'FontWeight','bold');
	xlabel(handles.axes3, 'Masked Image', 'FontSize', fontSize, 'Color', 'w', 'FontWeight','bold');
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
guidata(hObject, handles);


%=====================================================================================================================================================
% --- Executes when selected object is changed in panYAxis.
function panYAxis_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panYAxis
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Original_Image;
global binsToSuppress;
global g_MaxThreshold;
global g_MinThreshold;
try
	% 		plotLinearYAxis = handles.radCount.Value;
	% 		plotLinearYAxis = handles.radLogCount.Value;
	% Plot the histogram
	PlotHistogram(handles, Original_Image, binsToSuppress);
	% Place red bars at the threshold locations.
	PlaceThresholdBar(handles, g_MinThreshold, g_MaxThreshold);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from panYAxis_SelectionChangeFcn()


%====================================================================================================================================================================
% --- Executes on button press in chkScaleIntensity.
function chkScaleIntensity_Callback(hObject, eventdata, handles)
% hObject    handle to chkScaleIntensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of chkScaleIntensity
% Apply the threshold to the display so we can see its effect.
ShowThresholdedBinaryImage(hObject, eventdata, handles);

%====================================================================================================================================================================
% Pop up a warning message box with the warning message.
function WarnUser(warningMessage)
fprintf('%s\n', warningMessage);	% Also print it to the command window.
uiwait(warndlg(warningMessage));	% Pop up a message box
return; % from WarnUser()


%====================================================================================================================================================================
% --- Executes on button press in chkGridLines.
function chkGridLines_Callback(hObject, eventdata, handles)
global Original_Image;
global binsToSuppress;
global g_MaxThreshold;
global g_MinThreshold;
try
	% Plot the histogram
	PlotHistogram(handles, Original_Image, binsToSuppress);
	% Place red bars at the threshold locations.
	PlaceThresholdBar(handles, g_MinThreshold, g_MaxThreshold);
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from chkGridLines_Callback()


%====================================================================================================================================================================
function SetUpImpixelinfo(handles)
try
	hp = impixelinfo;
	hp.Units = 'normalized'; % Coordinates of the containing panel only, not the text inside.
% 	op = hp.Position;
	hp.Position = [0, 0.45, 0.2, 0.1];
	childObjects = findobj(hp);
	% childObjects(1) is the containing panel.
	% childObjects(2) is the label on the left that says "Pixel info:".
	% childObjects(3) is the string that says "(X, Y) Intensity" with the actual numbers.
	% Set up foreground color to be white.
	childObjects(2).ForegroundColor = [1,1,1];
	childObjects(3).ForegroundColor = [1,1,1];
	% Set up the font size.
	childObjects(2).FontSize = 11;
	childObjects(3).FontSize = 11;
	% Set up background color to be black.
	childObjects(2).BackgroundColor = [0,0,0]; % [R, G, B] in range 0-1
	childObjects(3).BackgroundColor = [0,0,0]; % [R, G, B] in range 0-1
	% Position the status label on the left, half way down.
	
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from SetUpImpixelinfo()


%=====================================================================================================================================================
% --- Executes on button press in btnTriangle.
function btnTriangle_Callback(hObject, eventdata, handles)
global g_MinThreshold;
global g_MaxThreshold;
try
	message = sprintf('Do you want bright objects or dark objects?');
	reply = questdlg(message, 'Toolbox missing', 'Bright', 'Dark', 'Bright');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
	theImage = getimage(handles.axes1);
	[counts, edges] = histcounts(theImage(:));
	if contains(reply, 'Dark')
		index = triangle_threshold(counts, 'L', false);
		thresholdValue = edges(index);
		handles.slider_Min.Value = handles.slider_Min.Min;
		handles.slider_Max.Value = thresholdValue;
		% Assign global variables so the threshold bars will be drawn at the correct locations.
		g_MinThreshold = handles.slider_Min.Value;
		g_MaxThreshold = thresholdValue;
	else
		index = triangle_threshold(counts, 'R', false);
		thresholdValue = edges(index);
		handles.slider_Min.Value = thresholdValue;
		handles.slider_Max.Value = handles.slider_Max.Max;
		% Assign global variables so the threshold bars will be drawn at the correct locations.
		g_MinThreshold = thresholdValue;
		g_MaxThreshold = handles.slider_Max.Value;
	end
	% Both thresholds have changed so we need to update both sliders and both edit fields.
	% Place the threshold bars and update the edit fields.
	slider_Max_Callback(hObject, eventdata, handles);
	slider_Min_Callback(hObject, eventdata, handles);
	
catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from btnTriangle_Callback()


%=====================================================================================================================================================
% --- Executes on button press in btnOtsu.
function btnOtsu_Callback(hObject, eventdata, handles)
global g_MinThreshold;
global g_MaxThreshold;
try
	message = sprintf('Do you want bright objects or dark objects?');
	reply = questdlg(message, 'Toolbox missing', 'Bright', 'Dark', 'Bright');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
	theImage = getimage(handles.axes1);
	thresholdValue = graythresh(theImage) * max(theImage(:));
	if contains(reply, 'Dark')
		handles.slider_Min.Value = handles.slider_Min.Min;
		handles.slider_Max.Value = thresholdValue;
		% Assign global variables so the threshold bars will be drawn at the correct locations.
		g_MinThreshold = thresholdValue;
		g_MaxThreshold = handles.slider_Max.Value;
	else
		handles.slider_Min.Value = thresholdValue;
		handles.slider_Max.Value = handles.slider_Max.Max;
		% Assign global variables so the threshold bars will be drawn at the correct locations.
		g_MinThreshold = thresholdValue;
		g_MaxThreshold = handles.slider_Max.Value;
	end
	% Both thresholds have changed so we need to update both sliders and both edit fields.
	% Place the threshold bars and update the edit fields.
	slider_Max_Callback(hObject, eventdata, handles);
	slider_Min_Callback(hObject, eventdata, handles);

catch ME
	errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
		ME.stack(1).name, ME.stack(1).line, ME.message);
	WarnUser(errorMessage);
end
return; % from btnOtsu_Callback()

	
