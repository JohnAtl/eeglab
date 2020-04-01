function varargout = PowPowCAT(varargin)
% POWPOWCAT MATLAB code for PowPowCAT.fig
%      POWPOWCAT, by itself, creates a new POWPOWCAT or raises the existing
%      singleton*.
%
%      H = POWPOWCAT returns the handle to a new POWPOWCAT or the handle to
%      the existing singleton*.
%
%      POWPOWCAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POWPOWCAT.M with the given input arguments.
%
%      POWPOWCAT('Property','Value',...) creates a new POWPOWCAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PowPowCAT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PowPowCAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PowPowCAT

% Last Modified by GUIDE v2.5 05-Dec-2017 17:03:44

% History
% 11/29/2017 Makoto. Incorporated Nattapoing's changes. Renamed the function to PowPowCAT for publication.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PowPowCAT_OpeningFcn, ...
                   'gui_OutputFcn',  @PowPowCAT_OutputFcn, ...
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


% --- Executes just before PowPowCAT is made visible.
function PowPowCAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PowPowCAT (see VARARGIN)

% Check data status
EEG = evalin('base', 'EEG');

if isempty(EEG.icaact)
    error(sprintf(...
        'EEG.icaact is empty. If you have not comptued ICA, run it first.\nAlso, check the checkbox for the option under ''File'' -> ''Memory and other options'' -> ''If set, precompute ICA activations...'''))
end

if size(EEG.icaact)==3
    error('Use continuous data.')
end

if isa(EEG.icaact, 'single')
    error(sprintf('EEG.icaact is single precision floating-point numeric array. Convert it to double precision.\nUncheck the checkbox ''If set, use single precision under Matlab 7.X...'' that can be found in ''File'' -> ''Memmory and other options'''))
end

% Pass the size of the covariance matrix to axes handle
if isfield(EEG.etc, 'PowPowCAT')
    freqData.freqs       = EEG.etc.PowPowCAT.freqs;
    freqData.numFreqBins = length(EEG.etc.PowPowCAT.freqs);
    set(hObject, 'UserData', freqData);
    disp('Pre-computed cross-frequency power spectum detected.')
else
    disp('Cross-frequency power spectum is not pre-computed yet.')
end

% Enter upper Frequency Limit
if isfield(EEG.etc, 'PowPowCAT')
    set(handles.upperFreqLimitEdit, 'String', [num2str(EEG.etc.PowPowCAT.freqs(end)) ' (calculated)']);
else
    set(handles.upperFreqLimitEdit, 'String', num2str(floor(EEG.srate/2)));
end

% Display no data logo on connectivityAxes.
logoPath = which('sccnLogo.jpg');
logoData = imread(logoPath);
logoData = mean(double(logoData),3);
clockOutputs = clock;
currentMonth = clockOutputs(2);
colorRatio = 0.12221946;
if currentMonth>=3 & currentMonth<=5
    originalColorMap = colormap(spring(256));
elseif currentMonth>=6 & currentMonth<=9
    originalColorMap = colormap(summer(256));
elseif currentMonth>=10 & currentMonth<=11
    originalColorMap = colormap(autumn(256)); 
else
    originalColorMap = colormap(winter(256));
end
customColorMap = originalColorMap*colorRatio + repmat(1-colorRatio, size(originalColorMap));
colormap(customColorMap)
axes(handles.covMatrixAxes)
imagesc(logoData);
set(gca, 'XTick', [], 'YTick', [])

% Display no data logo on specPlotAxes.
axes(handles.specPlotAxes)
imagesc(logoData);
set(gca, 'XTick', [], 'YTick', [])

% Display no data logo on corrPlotAxes.
set(handles.corrPlotAxes, 'Units', 'pixels');
corrPlotAxesSize = get(handles.corrPlotAxes, 'Position'); % Use 3rd and 4th for vertical and horizontal lengths
rectangleLogo1 = repmat(logoData(1,1), [length(logoData) round(length(logoData)*corrPlotAxesSize(3)/corrPlotAxesSize(4))]);
startPoint    = round(size(rectangleLogo1,2)/2)-round(size(logoData,2)/2);
rectangleLogo1(:,startPoint:startPoint+length(logoData)-1) = logoData;
axes(handles.corrPlotAxes)
imagesc(rectangleLogo1);
set(gca, 'XTick', [], 'YTick', [])

% Display no data logo on timeSeriesAxes.
set(handles.timeSeriesAxes, 'Units', 'pixels');
timeSeriesAxesSize = get(handles.timeSeriesAxes, 'Position'); % Use 3rd and 4th for vertical and horizontal lengths
rectangleLogo2 = repmat(logoData(1,1), [length(logoData) round(length(logoData)*timeSeriesAxesSize(3)/timeSeriesAxesSize(4))]);
startPoint    = round(size(rectangleLogo2,2)/2)-round(size(logoData,2)/2);
rectangleLogo2(:,startPoint:startPoint+length(logoData)-1) = logoData;
axes(handles.timeSeriesAxes)
imagesc(rectangleLogo2);
set(gca, 'XTick', [], 'YTick', [])

% Display no data logo on timeSeriesScatterAxes.
axes(handles.timeSeriesScatterAxes)
imagesc(rectangleLogo2);
set(gca, 'XTick', [], 'YTick', [])

drawnow

% Change the figure title
set(gcf, 'Name', 'PowPowCAT()');

% Choose default command line output for PowPowCAT
handles.output = hObject;

% % Show toolbar
% set(gcf, 'ToolBar', 'figure')

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = PowPowCAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function icIdxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to icIdxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of icIdxEdit as text
%        str2double(get(hObject,'String')) returns contents of icIdxEdit as a double

EEG = evalin('base', 'EEG');

% Plot spectrum with SD
currentIcIdx = str2num(get(handles.icIdxEdit, 'String'));
currentSpec  = 10*log10(EEG.etc.PowPowCAT.meanPSD(:,currentIcIdx));
freqs = EEG.etc.PowPowCAT.freqs;
axes(handles.specPlotAxes)
cla
plot(currentSpec, 'LineWidth', 2)

% covMatix = evalin('base', 'EEG.etc.PowPowCAT.covMatrix');
userData.currentIcIdx= currentIcIdx;
userData.currentSpec = currentSpec;
userData.freqs       = freqs;
userData.corrCoeff   = EEG.etc.PowPowCAT.covMatrix(:,:, currentIcIdx);
userData.minValue    = min(vec(userData.corrCoeff));
% Natty ----------- Adding
userData.time          = EEG.xmax/200:EEG.xmax/200:EEG.xmax; % sec
userData.timeSeriesPSD = EEG.etc.PowPowCAT.timeSeriesPSD(:,:,currentIcIdx);
% -----------------
set(handles.specPlotAxes, 'UserData', userData);

% http://www.ibm.com/support/knowledgecenter/ja/SSLVMB_22.0.0/com.ibm.spss.statistics.help/vis_workbench/graphboard_editing_axes.htm
% safelog == sign(x)*log(1+abs(x))
if get(handles.plotSdCheckbox, 'value')==1
    hold on
    upperEnv     = 10*log10(EEG.etc.PowPowCAT.meanPSD(:,currentIcIdx) + EEG.etc.PowPowCAT.stdPSD(:,currentIcIdx));
    lowerEnv     = currentSpec-(upperEnv-currentSpec);
    fillHandle = fill([1:length(freqs) length(freqs):-1:1], [upperEnv' lowerEnv(end:-1:1)'], [0.66 0.76 1], 'Linestyle', 'none');
    uistack(fillHandle, 'bottom');
    hold off
end
XTickLabels = round(freqs(10:10:length(freqs))*10)/10;
XTickLabels(XTickLabels>10) = round(XTickLabels(XTickLabels>10));
set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', XTickLabels, 'fontsize', 11);
% Natty ----------- Disabling
%xlabel('Frequency (Hz)', 'fontsize', 12)
% -----------------
ylabel(sprintf('Power Spectrum Density by spectrogram()\n10*log10(uV^2)/Hz'), 'fontsize', 12)



% Prepare covariance matrix and mask mask if requested.
currentIcCovMatrix = EEG.etc.PowPowCAT.covMatrix(:,:,currentIcIdx);
if     get(handles.plotTypePopupmenu, 'value')==1
    currentMask = ones(size(currentIcCovMatrix));
elseif get(handles.plotTypePopupmenu, 'value')==2
    currentMask = EEG.etc.PowPowCAT.pvalMatrix(:,:,currentIcIdx);
    currentMask = currentMask<EEG.etc.PowPowCAT.fdr005;
elseif get(handles.plotTypePopupmenu, 'value')==3
    currentMask = EEG.etc.PowPowCAT.pvalMatrix(:,:,currentIcIdx);
    currentMask = currentMask<EEG.etc.PowPowCAT.fdr001;
end



% Plot correlation matrix
axes(handles.covMatrixAxes)
imagesc(currentIcCovMatrix.*currentMask, [-0.8 0.8]);
customColorMap = colormap(jet);
if get(handles.plotTypePopupmenu, 'value')~=1
    customColorMap(129,:) = [0.5 0.5 0.5]; % Mask exact zeros with gray pixels
end
colormap(customColorMap)
currentAxesPosition = get(gca, 'position');
colorbarHandle = colorbar;
set(get(colorbarHandle, 'title'), 'String', 'Corr. Coef')
set(handles.covMatrixAxes, 'position', currentAxesPosition);
axis xy
axis square
tickLabels = round(freqs(10:10:length(freqs))*10)/10;
tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', tickLabels,...
         'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels,...
         'fontsize', 12)
xlabel('Frequency (Hz)', 'fontsize', 12)
ylabel('Frequency (Hz)', 'fontsize', 12)
% Natty ----------- Spectrum --> Spectral
title('Spectral Covariance (Click!)', 'fontsize', 18, 'Color', [0 0 0.4])
% -----------------

% Enable the frequency range report
set(handles.freqCoupleDisplayText, 'enable', 'on')

% Choose default command line output for PowPowCAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



    






% --- Executes during object creation, after setting all properties.
function icIdxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to icIdxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in precomputePushbutton.
function precomputePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to precomputePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check data status
EEG = evalin('base', 'EEG');

%     Spectrogram using a Short-Time Fourier Transform (STFT).
%     S = spectrogram(X) returns the spectrogram of the signal specified by
%     vector X in the matrix S. By default, X is divided into eight segments
%     with 50% overlap, each segment is windowed with a Hamming window. The
%     number of frequency points used to calculate the discrete Fourier
%     transforms is equal to the maximum of 256 or the next power of two
%     greater than the length of each segment of X.
%
%     Each column of S contains an estimate of the short-term, time-localized
%     frequency content of the signal X.  Time increases across the columns
%     of S, from left to right.  Frequency increases down the rows, starting
%     at 0.  If X is a length NX complex signal, S is a complex matrix with
%     NFFT rows and k = fix((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP)) columns.
%     For real X, S has (NFFT/2+1) rows if NFFT is even, and (NFFT+1)/2 rows
%     if NFFT is odd.  
%  
%     [S,F,T] = spectrogram(...) returns a vector of frequencies F and a
%     vector of times T at which the spectrogram is computed. F has length
%     equal to the number of rows of S. T has length k (defined above) and
%     its value corresponds to the center of each segment.
%  
%     [S,F,T] = spectrogram(X,WINDOW,NOVERLAP,F,Fs) where F is a vector of 
%     frequencies in Hz (with 2 or more elements) computes the spectrogram at 
%     those frequencies using the Goertzel algorithm. The specified 
%     frequencies in F are rounded to the nearest DFT bin commensurate with 
%     the signal's resolution. 
%
%     [S,F,T,P] = spectrogram(...) P is a matrix representing the Power
%     Spectral Density (PSD) of each segment. For real signals, spectrogram
%     returns the one-sided modified periodogram estimate of the PSD of each
%     segment; for complex signals and in the case when a vector of
%     frequencies is specified, it returns the two-sided PSD.  

% Compute PSD with 1-sec window with 50% overlap
% Natty ----------- Disabling and Changing to more linear log transform.
%freqBins = logspace(0, log10(str2num(get(handles.upperFreqLimitEdit, 'String'))), 100);
deviationFromLog = 5;
freqBins = logspace(log10(1+deviationFromLog), log10(str2num(get(handles.upperFreqLimitEdit, 'String'))+deviationFromLog), 100)-deviationFromLog; %5th
% -----------------
processTimeList = zeros(size(EEG.icaact,1),1);
for icIdx = 1:size(EEG.icaact,1)
    
    % Compute short-term Fourier Transform
    % Obtain the PSD length from the output of the first iteration (Thanks Pal Gunnar Lasson! 1/10/2017)
    tic;
    if icIdx == 1;
        [~, freqs, times, firstPSD] = spectrogram(EEG.icaact(icIdx,:), EEG.srate, floor(EEG.srate/2), freqBins, EEG.srate);
        PSD = zeros(100, size(firstPSD,2), size(EEG.icaact,1));
        PSD(:,:,icIdx) = firstPSD;
    else
        [~, ~, ~, PSD(:,:,icIdx)] = spectrogram(EEG.icaact(icIdx,:), EEG.srate, floor(EEG.srate/2), freqBins, EEG.srate);
    end
    
    % Display time elapsed
    timeElapsed = toc;
    processTimeList(icIdx) = timeElapsed;
    if icIdx >= 4
        robustIdx = find(processTimeList(1:icIdx)>=prctile(processTimeList(1:icIdx),25) & processTimeList(1:icIdx)<=prctile(processTimeList(1:icIdx),75));
        meanValueSoFar = mean(processTimeList(robustIdx));
    else
        meanValueSoFar = mean(processTimeList(1:icIdx));
    end
    timeRemaining = (length(processTimeList)-icIdx)*meanValueSoFar;
    if     timeRemaining >= 3600
        disp(sprintf('%.0f/%.0f done, %.1f hours remaining.',   icIdx, size(EEG.icaact,1), timeRemaining/3600));
    elseif timeRemaining >  60
        disp(sprintf('%.0f/%.0f done, %.1f minutes remaining.', icIdx, size(EEG.icaact,1), timeRemaining/60));
    else
        disp(sprintf('%.0f/%.0f done, %.0f seconds remaining.', icIdx, size(EEG.icaact,1), timeRemaining));
    end
end



% Remove the chunks that contain 'boundary'
if isempty(EEG.event)
    boundaryIdx = [];
else
    boundaryIdx = find(strcmp({EEG.event.type}, 'boundary'));
end
if any(boundaryIdx)
    boundaryLatencyInSec = [EEG.event(boundaryIdx).latency]*1/EEG.srate;
    [~, removeIdx] = histc(boundaryLatencyInSec, times-0.5); % 'times' are bin centers, so converting to bin edges.
else
    removeIdx = [];
end
goodTimesIdx = setdiff(1:length(times), removeIdx);
times = times(goodTimesIdx);
PSD   = PSD(:,goodTimesIdx,:);

% Natty ----------- Disabling
% Clean data.
% residualVariance = (bsxfun(@minus, PSD, median(PSD,2))).^2;
% residualVarianceSum = squeeze(sum(sum(residualVariance,3),1));
% goodChunkIdx = find(residualVarianceSum<=prctile(residualVarianceSum, 90)); % Discard 10% of chunks
% cleanPSD = PSD(:,goodChunkIdx,:);
%     % figure
%     % plot(residualVarianceSum)
%     % figure
%     % plot(residualVarianceSum(goodChunkIdx))
% -----------------
% Natty ----------- Adding
% Clean data by SD.
residualVariance2 = bsxfun(@minus, PSD, mean(PSD,2));
residualVariance_norm = bsxfun(@rdivide, residualVariance2, std(PSD,0,2));
residualVariance_norm_abs = abs(residualVariance_norm);
residualVarianceSum = squeeze(sum(sum(residualVariance_norm_abs,3),1));
goodChunkIdx2 = find(residualVarianceSum<=prctile(residualVarianceSum, 80));
cleanPSD = PSD(:,goodChunkIdx2,:);
clear PSD residualVariance
% -----------------

% Compute cross-frequency coupling for each ICs (spectrum covariance across moving windows).
covMatrix = zeros(size(cleanPSD,1), size(cleanPSD,1), size(cleanPSD,3));
for icIdx = 1:size(covMatrix,3)
    covMatrix(:,:,icIdx) = corrcoef(cleanPSD(:,:,icIdx)');
end

% Natty ----------- Adding
% Subsample the PSD time series into 200 points
cleanPsdPerm   = permute(cleanPSD, [2 1 3]);
cleanPsdPerm2D = cleanPsdPerm(:,:);
resampledPSD   = resample(cleanPsdPerm2D, 200, size(cleanPsdPerm2D,1));
resampledPsd3D = permute(reshape(resampledPSD, [200, 100, size(cleanPSD,3)]), [2 1 3]);
psdTimeSeriesNorm   = bsxfun(@rdivide, resampledPsd3D, std( resampledPsd3D,0,2));
psdTimeSeriesDemean = bsxfun(@minus,   psdTimeSeriesNorm, mean(psdTimeSeriesNorm,2));
clear cleanPsdPerm cleanPsdPerm2D resampledPSD resampledPsd3D psdTimeSeriesNorm
% -----------------

% Perform permutation test.
numIterations = 1000 %5000;
processTimeList = zeros(numIterations,1);
surroMatrix = zeros(size(cleanPSD,1), size(cleanPSD,1), size(cleanPSD,3), numIterations);
for iterationIdx = 1:numIterations
    tic;
    permIdxMatrix = zeros(size(cleanPSD,1), size(cleanPSD,2));
    for freqIdx = 1:size(cleanPSD,1)
        permIdxMatrix(freqIdx,:) = randperm(size(cleanPSD,2));
    end
    permIdxMatrix = repmat(permIdxMatrix, [1 1 size(cleanPSD,3)]);
    surroPSD = cleanPSD(permIdxMatrix);
    for icIdx = 1:size(surroMatrix,3)
        surroMatrix(:,:,icIdx,iterationIdx) = corrcoef(surroPSD(:,:,icIdx)');
    end
    
    % Display time elapsed.
    timeElapsed = toc;
    processTimeList(iterationIdx) = timeElapsed;
    if iterationIdx >= 4
        robustIdx = find(processTimeList(1:iterationIdx)>=prctile(processTimeList(1:iterationIdx),25) & processTimeList(1:iterationIdx)<=prctile(processTimeList(1:iterationIdx),75));
        meanValueSoFar = mean(processTimeList(robustIdx));
    else
        meanValueSoFar = mean(processTimeList(1:iterationIdx));
    end
    timeRemaining = (length(processTimeList)-iterationIdx)*meanValueSoFar;
    if mod(iterationIdx, 100) == 0
        if     timeRemaining >= 3600
            disp(sprintf('%.0f/%.0f done, %.1f hours remaining.',   iterationIdx, numIterations, timeRemaining/3600));
        elseif timeRemaining >  60
            disp(sprintf('%.0f/%.0f done, %.1f minutes remaining.', iterationIdx, numIterations, timeRemaining/60));
        else
            disp(sprintf('%.0f/%.0f done, %.0f seconds remaining.', iterationIdx, numIterations, timeRemaining));
        end
    end
end

% Perform Tim's non-parametric statistics.
pvalMatrix = zeros(size(covMatrix));
for icIdx = 1:size(covMatrix,3)
    disp(sprintf('Non-parametric test for each IC: %2.0f/%2.0f...', icIdx, size(pvalMatrix,3)))
    tmpSurro = squeeze(surroMatrix(:,:,icIdx,:));
    tmpCov   = covMatrix(:,:,icIdx);
    pvalMatrix(:,:,icIdx) = stat_surrogate_pvals(tmpSurro, tmpCov, 'both');
end

% FDR correction using all pixels from all ICs for global correction (this one is a very easy stats...)
pvalMask = pvalMatrix;
pvalMatrix_noDiagonal = pvalMatrix;
pvalMatrix_noDiagonal = bsxfun(@plus, pvalMatrix_noDiagonal, 2*eye(length(freqBins)));
nonDiagonalPvalues    = pvalMatrix_noDiagonal(pvalMatrix_noDiagonal<=1); % Must be integer multiple of 100^2-100
fdr005 = fdr(nonDiagonalPvalues, 0.05);
fdr001 = fdr(nonDiagonalPvalues, 0.01);

% Store the result.
EEG.etc.PowPowCAT.covMatrix  = covMatrix;
EEG.etc.PowPowCAT.freqs      = freqs;
EEG.etc.PowPowCAT.meanPSD    = squeeze(mean(cleanPSD,2));  % dB-converted
EEG.etc.PowPowCAT.stdPSD     = squeeze(std(cleanPSD,0,2)); % dB-converted
EEG.etc.PowPowCAT.pvalMatrix = pvalMatrix;
EEG.etc.PowPowCAT.fdr005 = fdr005;
EEG.etc.PowPowCAT.fdr001 = fdr001;
% Natty -----------
EEG.etc.PowPowCAT.timeSeriesPSD = psdTimeSeriesDemean;
% -----------------
% Store the resutls to handles.PowPowCATFigure, 'UserData'
freqData.freqs       = EEG.etc.PowPowCAT.freqs;
freqData.numFreqBins = length(EEG.etc.PowPowCAT.freqs);
set(handles.PowPowCATFigure, 'UserData', freqData);

% Save EEG.
assignin('base', 'EEG', EEG)

% Choose default command line output for PowPowCAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% Natty -----------
disp('Pre-computing: Done');
% -----------------




% --- Executes on selection change in plotTypePopupmenu.
function plotTypePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotTypePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotTypePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotTypePopupmenu

icIdxEdit_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function plotTypePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotTypePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperFreqLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to upperFreqLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperFreqLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of upperFreqLimitEdit as a double


% --- Executes during object creation, after setting all properties.
function upperFreqLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperFreqLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Plot35IcPushbutton.
function Plot35IcPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Plot35IcPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EEG = evalin('base', 'EEG');
freqs = EEG.etc.PowPowCAT.freqs;

figure
set(gcf, 'color', [0.66 0.76 1])
for icIdx = 1:35
    
    % Plot correlation matrix
    subplot(5,7,icIdx)
    imagesc(EEG.etc.PowPowCAT.covMatrix(:,:,icIdx), [-0.8 0.8]);
    colormap(jet);
    axis xy
    axis square
    
    % Zoom in a little bit
    gcaPosition = get(gca, 'position');
    gcaPosition([3 4]) = gcaPosition([3 4])*1.1;
    set(gca, 'position', gcaPosition);
    tickLabels = round(freqs(10:10:length(freqs))*10)/10;
    tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
    set(gca, 'XTick', [], 'XTickLabel', [],...
        'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels,...
        'fontsize', 10)
    title(['IC ' num2str(icIdx)], 'fontsize', 14)
    if icIdx == size(EEG.icaweights,1)
        return
    end
end


% --- Executes on mouse motion over figure - except title and menu.
function PowPowCATFigure_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to PowPowCATFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do not activate until pre-compute is done.
if isempty(get(handles.PowPowCATFigure, 'UserData'))
    return
end

% Obtain cursor position within the Axes.
cursorCoordinate = get(handles.covMatrixAxes, 'CurrentPoint');
        % % display coordinate
        % disp(sprintf('row==%d column=%d', round(cursorCoordinate(1,2)), round(cursorCoordinate(1,1))));

freqData = get(handles.PowPowCATFigure, 'UserData');
numFreqBins = freqData.numFreqBins;
% Show the current combination of freqs
if cursorCoordinate(1,1)>0.5 & cursorCoordinate(1,1)<numFreqBins+0.5 & cursorCoordinate(1,2)>0.5 & cursorCoordinate(1,2)<numFreqBins+0.5

    %PowPowCAT = evalin('base', 'EEG.etc.PowPowCAT');
    
    % Obtain freq index.
    freq1Idx = round(cursorCoordinate(1,1));
    freq2Idx = round(cursorCoordinate(1,2));
    
    % Load current IC data attached to the PSD axes
    currentIcData = get(handles.specPlotAxes, 'UserData');
    if isempty(currentIcData)
        return
    end
    
    % Display the coupled frequencies.
    corrCoeff = currentIcData.corrCoeff(freq2Idx, freq1Idx);
    freqs     = currentIcData.freqs;
    str       = sprintf('%.1fHz with %.1fHz, r=%.3f', freqs(freq2Idx), freqs(freq1Idx), corrCoeff);
    set(handles.freqCoupleDisplayText, 'String', str, 'fontsize', 16)

    % Plot one-to-multi-freq correlation coeff plot
    axes(handles.corrPlotAxes)
    plot(currentIcData.corrCoeff(freq2Idx,:), 'linewidth', 2)
    XTickLabels = round(freqs(10:10:length(freqs))*10)/10;
    XTickLabels(XTickLabels>10) = round(XTickLabels(XTickLabels>10));
    set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', XTickLabels, 'fontsize', 11);
    line([freq2Idx freq2Idx], ylim, 'color', [1 0 0])
    % Natty ----------- Disabling
    %ylim([currentIcData.minValue 0.8])
    % -----------------
    % Natty ----------- Adding
    ylim([0 1])
    set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]) 
    xlabel('Frequency (Hz)')
    % -----------------
    ylabel('Corr. Coeff. (r)')
    % Natty ----------- Adding
    drawnow
    % -----------------
    
%     % Plot vertical bars on spectraplot
%     axes(handles.specPlotAxes);
%     axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
%     delete(axesHandlesToChildObjects);
%     hold on
%     plot(currentIcData.currentSpec, 'LineWidth', 2)
%     line([freq2Idx freq2Idx], ylim, 'color', [1 0 0]);
%     line([freq1Idx freq1Idx], ylim);
%     hold off
end


% --- Executes on button press in plotSdCheckbox.
function plotSdCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to plotSdCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotSdCheckbox

icIdxEdit_Callback(hObject, eventdata, handles)



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function PowPowCATFigure_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to PowPowCATFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do not activate until pre-compute is done.
if isempty(get(handles.PowPowCATFigure, 'UserData'))
    return
end

% Obtain cursor position within the Axes.
cursorCoordinate = get(handles.covMatrixAxes, 'CurrentPoint');
        % % display coordinate
        % disp(sprintf('row==%d column=%d', round(cursorCoordinate(1,2)), round(cursorCoordinate(1,1))));

freqData = get(handles.PowPowCATFigure, 'UserData');
numFreqBins = freqData.numFreqBins;
% Show the current combination of freqs
if cursorCoordinate(1,1)>0.5 & cursorCoordinate(1,1)<numFreqBins+0.5 & cursorCoordinate(1,2)>0.5 & cursorCoordinate(1,2)<numFreqBins+0.5

    %PowPowCAT = evalin('base', 'EEG.etc.PowPowCAT');
    
    % Obtain freq index.
    freq1Idx = round(cursorCoordinate(1,1));
    freq2Idx = round(cursorCoordinate(1,2));
    
    % Load current IC data attached to the PSD axes
    currentIcData = get(handles.specPlotAxes, 'UserData');
    if isempty(currentIcData)
        return
    end
    
    % Display the coupled frequencies.
    corrCoeff = currentIcData.corrCoeff(freq2Idx, freq1Idx);
    freqs     = currentIcData.freqs;
    str       = sprintf('%.1fHz with %.1fHz, r=%.3f', freqs(freq2Idx), freqs(freq1Idx), corrCoeff);
    set(handles.freqCoupleDisplayText, 'String', str, 'fontsize', 16)

    % Plot one-to-multi-freq correlation coeff plot
    axes(handles.corrPlotAxes)
    plot(currentIcData.corrCoeff(freq2Idx,:), 'linewidth', 2)
    XTickLabels = round(freqs(10:10:length(freqs))*10)/10;
    XTickLabels(XTickLabels>10) = round(XTickLabels(XTickLabels>10));
    set(gca, 'XTick', 10:10:length(freqs), 'XTickLabel', XTickLabels, 'fontsize', 11);
    line([freq2Idx freq2Idx], ylim, 'color', [1 0 0])
    % Natty ----------- Disabling
    %ylim([currentIcData.minValue 0.8])
    % -----------------
    % Natty ----------- Adding
    ylim([0 1])
    set(gca, 'YTick', [0 0.2 0.4 0.6 0.8]) 
    xlabel('Frequency (Hz)')
    % -----------------
    ylabel('Corr. Coeff. (r)')

    % Plot vertical bars on spectraplot
    axes(handles.specPlotAxes);
    axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
    delete(axesHandlesToChildObjects);
    hold on
    plot(currentIcData.currentSpec, 'LineWidth', 2)
    line([freq2Idx freq2Idx], ylim, 'color', [1 0 0]);
    line([freq1Idx freq1Idx], ylim);
    hold off
    
    % Natty ----------- Adding
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot two PSD time-series. %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(handles.timeSeriesAxes)
    cla
    timeSeries1 = currentIcData.timeSeriesPSD(freq1Idx,:);
    timeSeries2 = currentIcData.timeSeriesPSD(freq2Idx,:);
    plot(currentIcData.time, timeSeries2, 'r', 'linewidth', 1); hold on
    plot(currentIcData.time, timeSeries1, 'linewidth', 1);
    xlim([currentIcData.time(1) currentIcData.time(end)])
    ylim([min([timeSeries1 timeSeries2]) max([timeSeries1 timeSeries2])])
    xlabel('Time (s)', 'fontsize', 11)
    ylabel('Normalized power','fontsize', 11)
    
    % Plot PSD time-series
    axes(handles.timeSeriesScatterAxes)
    cla
    corrCoeff = currentIcData.corrCoeff(freq2Idx, freq1Idx);
    scatter(timeSeries1,timeSeries2, 10, [0.66 0.76 1]); hold on
    line([mean(timeSeries1) mean(timeSeries1)],[min(timeSeries2) max(timeSeries2)],'Color','k','LineStyle','-');
    line([min(timeSeries1) max(timeSeries1)],[mean(timeSeries2) mean(timeSeries2)],'Color','k','LineStyle','-');
    line([min(timeSeries1) max(timeSeries1)],[min(timeSeries1)*corrCoeff max(timeSeries1)*corrCoeff],'Color',[1,0,0],'LineStyle','-');
    xlim([min(timeSeries1) max(timeSeries1)]);
    ylim([min(timeSeries2) max(timeSeries2)]);
    xlabel('Selected freq. pow. 1', 'fontsize', 11)
    ylabel('Selected freq. pow. 2', 'fontsize', 11)
    %title(['r = ' num2str(corr_val)],'fontsize', 10);
    hold off
    % -----------------
end


% --- Executes when PowPowCATFigure is resized.
function PowPowCATFigure_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to PowPowCATFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
