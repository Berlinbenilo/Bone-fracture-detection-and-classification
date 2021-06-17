function varargout = fracture_detection(varargin)
% FRACTURE_DETECTION MATLAB code for fracture_detection.fig
%      FRACTURE_DETECTION, by itself, creates a new FRACTURE_DETECTION or raises the existing
%      singleton*.
%
%      H = FRACTURE_DETECTION returns the handle to a new FRACTURE_DETECTION or the handle to
%      the existing singleton*.
%
%      FRACTURE_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRACTURE_DETECTION.M with the given input arguments.
%
%      FRACTURE_DETECTION('Property','Value',...) creates a new FRACTURE_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fracture_detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fracture_detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fracture_detection

% Last Modified by GUIDE v2.5 17-Jun-2021 07:37:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fracture_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @fracture_detection_OutputFcn, ...
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


% --- Executes just before fracture_detection is made visible.
function fracture_detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fracture_detection (see VARARGIN)

% Choose default command line output for fracture_detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

han = axes('unit','normalized','position',[0 0 1 1]);
new = imread('A.jpg');
imagesc(new);
set(han,'handlevisibility','on','visible','off');
uistack(han,'top');

% UIWAIT makes fracture_detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fracture_detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global img imag
load detector1.mat
[filename,filepath] = uigetfile('file selector');
fullpath = strcat(filepath,filename);
imag = imread(fullpath);
img = imag;
axes(handles.axes1);
imshow(img);
title('Input image');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global img boneEdges
ImgBlurSigma = 2; % Amount to denoise input image

img = (rgb2gray(img)); % Load image
img = imfilter(img, fspecial('gaussian', 10, ImgBlurSigma), 'symmetric'); % Denoise
boneEdges = edge(img, 'canny');
axes(handles.axes2)
imshow(boneEdges)
title('Edge detection canny');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global boneEdges img
MinHoughPeakDistance = 5; % Distance between peaks in Hough transform angle detection
HoughConvolutionLength = 40; % Length of line to use to detect bone regions
HoughConvolutionDilate = 2; % Amount to dilate kernel for bone detection
BreakLineTolerance = 0.25; % Tolerance for bone end detection
breakPointDilate = 6; % Amount to dilate detected bone end points
boneEdges = bwmorph(boneEdges, 'close');
edgeRegs = regionprops(boneEdges, 'Area', 'PixelIdxList');
AreaList = sort(vertcat(edgeRegs.Area), 'descend');
edgeRegs(~ismember(vertcat(edgeRegs.Area), AreaList(1:2))) = [];
edgeImg = zeros(size(img, 1), size(img,2));
edgeImg(vertcat(edgeRegs.PixelIdxList)) = 1;
[H,T,R] = hough(edgeImg,'RhoResolution',1,'Theta',-90:2:89.5);
maxHough = max(H, [], 1);
HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);
[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);

% Plot Hough detection results
axes(handles.axes2)
plot(T, maxHough);
hold on
plot([min(T) max(T)], [HoughThresh, HoughThresh], 'r');
plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold off
xlabel('Theta Value'); ylabel('Max Hough Transform');
legend({'Max Hough Transform', 'Hough Peak Threshold', 'Detected Peak'});

% Locate site of break
if numel(HoughPeaks) > 1;
    BreakStack = zeros(size(img, 1), size(img, 2), numel(HoughPeaks));
    % Convolute edge image with line of detected angle from hough transform
    for m = 1:numel(HoughPeaks);

        boneKernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));
        kern = double(bwmorph(boneKernel.getnhood(), 'dilate', HoughConvolutionDilate));
        BreakStack(:,:,m) = imfilter(edgeImg, kern).*edgeImg;
    end

    % Take difference between convolution images.  Where this crosses zero
    % (within tolerance) should be where the break is.  Have to filter out
    % regions elsewhere where the bone simply ends.
    brImg = abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0;
    [BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0);
    brImg = bwmorph(brImg, 'dilate', breakPointDilate);
    brReg = regionprops(brImg, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
        'Orientation', 'Centroid');
    brReg(vertcat(brReg.Area) ~= max(vertcat(brReg.Area))) = [];

    % Calculate bounding ellipse
    brReg.EllipseCoords = zeros(100, 2);
    t = linspace(0, 2*pi, 100);
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    brReg = [];

end

% Draw ellipse around break location

% figure(2)
% imshow(img)
% hold on
% colormap('gray')
% if ~isempty(brReg)
%     plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
% end
% hold off



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global img
% digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
%     'nndatasets','DigitDataset');
% imds = imageDatastore(digitDatasetPath, ...
%     'IncludeSubfolders',true,'LabelSource','foldernames');

dataset = 'C:/Users/BERLIN/Desktop/Bone fracture detection/Image';
imds = imageDatastore(dataset, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

labelCount = countEachLabel(imds);
img = readimage(imds,1);
size(img)
numTrainFiles = 80;
[imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');

layers = [
    imageInputLayer([300 300 3])
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',10, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');

net = trainNetwork(imdsTrain,layers,options);
% save Mynetwork.mat net
% 
% 
% % image = imread(imdsValidation.Files{1,1});
% load Mynetwork.mat
image = img;
out=classify(net,image);
c = cellstr(out)
msgbox(sprintf(c{1,1}))


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global imag

load detector1.mat
[bboxes,scores] = detect(acfDetector,imag);

flag=0;
for i = 1:length(scores)
    if scores(i,1) >= 85.2 && scores(i,1)<= 91
        flag=1;
        annotation = sprintf('Fracture: %.1f',scores(i));
        imag = insertObjectAnnotation(imag,'rectangle',bboxes(i,:),annotation,'LineWidth',4);
%         flag = flag+1;
%         bbox1(flag,1:4) = bboxes(i,1:4);
%         score1(flag,1) = scores(i,1);
%         boundingBoxArea(flag,1) = bbox1(flag,3) * bbox1(flag,4);
    end
end
axes(handles.axes2)
imshow(imag);
title('Detection area')
if flag==0
    msgbox('No fracture detected')
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% arrayfun(@cla,findall(0,'type','axes'))
set(handles.axes1,'visible','off')
set(handles.axes2,'visible','off')
cla(handles.axes1)
cla(handles.axes2)