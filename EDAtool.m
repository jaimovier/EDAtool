%% EDA Signal Tool
% Javier Jaimovich (2012)
% version 1.6 (updated 15-08-2012)
% 
% [EDA_f, EDA_x, EDA_P, EDA_T, Q] = ...
%     EDAtool(EDA, SR, debug, InterOption, ResampleOption, par)
%
% Returns low-pass version of the signal & removes artefacts, with the
% option of interpolating between the artefact. Phasic (fast changes) and 
% Tonic (slow changes) components are obtained via IIR filtering. The 
% function resamples the EDA signal to 50Hz, with the option of 
% resampling back to the original SR.
%
% EDA_f: low-pass filtered EDA
% EDA_x: EDA_f without artefacts (NaN or Interpolation)
% EDA_P: Phasic EDA or EDR (with NaN)
% EDA_T: Tonic EDA or EDL
% Q: Confidence of EDA signal (%)
%
% EDA: EDA signal
% SR: Signal's Sample Rate
% debug: 0 for no debugging, 1 to debug (includes plot)
% ResampleOption: 0 for no resampling, 1 to resample
% InterOption: Option to interpolate between artefacts (0:off, 1:on)
% par is a vector with [RangeMin RangeMax WindowSize Threshold NaNSize]
% being:
%   RangeMin: Lowest value in EDA range (e.g. open circuit)
%   RangeMax: Highest value in EDA range (e.g. closed circuit)
%   WindowSize: Size of artefact window in seconds.
%   Threshold: Threshold for detecting artefacts within window (0-100)
%   NaNSize: Size of NaN replacement vector for each artefact.
%   
%   Default values:
%       debug = 0
%       ResampleOption = 0
%       RangeMin = 150 (for BioControl sensor with 10bit ADC)
%       RangeMax = 490 (for BioControl sensor)
%       WindowSize = 0.2 [s]
%       Threshold = 4%
%       NaNSize = 1.5 [s]
%
% Note: Values outside or equal to range limits will be considered
% artefacts.

function [EDA_f, EDA_x, EDA_P, EDA_T, Q] = ...
    EDAtool(EDA, SR, debug, InterOption, ResampleOption, par)

if (debug==1);fprintf('EDA tool v1.6\n');end

%% DEFINITIONS (and default values)
SR_original = SR;
SR = 50; %Target sample rate for resampling and for saved filter coef.
art_buffer = 0.2; % buffer size for artefact detection (in seconds)
NaN_buffer = 1.5; %NaN window size (in seconds) for artefact replacement
EDA_range = [150 490]; %measured for BioControl circuit with 10bit ADC
thres = 4; %Default threshold for artefact detection

%% Check input parameters
if nargin < 2 || nargin >  6
    error('Check input parameters')
end
if SR <= 0
    error('SR must be greater than 0')
end
if length(EDA) < 5*SR
    error('EDA signal must be longer than 5 sec')
end
if exist('debug','var') == 0; debug = 0; end
if exist('ResampleOption','var') == 0; ResampleOption = 0; end
if exist('InterOption','var') == 0; InterOption = 1; end
if exist('par','var') == 1;
    if length(par)~=5; error('Par vector needs 5 parameters'); end
    EDA_range = par(1:2);
    art_buffer = par(3);
    thres = par(4);
    NaN_buffer = par(5);
    if (thres < 0)||(art_buffer<=0)||(EDA_range(1)>=EDA_range(2))||(NaN_buffer<=0)
        error('Check Parameters!'); end
end

if min(size(EDA))>1; error('EDA must be a vector Nx1 or 1xN');end

%% Pre-processing

N_o = length(EDA); %save original size for final resampling

%Resample to target SR (50Hz)
[P, Q] = rat(SR/SR_original); %ratio between SRs
EDA = resample(EDA-EDA(1),P,Q)+EDA(1);

%Logical vector with index of samples outside range
OutRange = logical((EDA<=EDA_range(1))+(EDA>=EDA_range(2)));

EDA = (EDA-EDA_range(1))/range(EDA_range); %Normalize
EDA = max(0,EDA); EDA = min(1,EDA); %limit to values within 0 1

EDA_o = EDA(1); % save 1st value
EDA = EDA - EDA_o; %Center EDA signal (for filters)

N = length(EDA); %length of resampled vector

%% Low Pass filter - Remove noise

%check if target SR has changed
if SR~=50; error('Target SR != from filter coefficients');end

%load filter coefficients
load('EDAtool_coeff_LP_50HzSR.mat');

% Response: 'Lowpass'
% Specification: 'Fp,Fst,Ap,Ast'
% Description: {4x1 cell}
% NormalizedFrequency: false
% Fs: 50
% Fpass: 0.5
% Fstop: 1
% Apass: 1
% Astop: 60
% 
% FilterStructure: 'Direct-Form II, Second-Order Sections'
% Arithmetic: 'double'
% sosMatrix: [1x6 double]
% ScaleValues: [3.94713593029321e-06;1]
% OptimizeScaleValues: true
% PersistentMemory: false

EDA_f = filter(Hlp_noise,EDA);

%Compensate filter latency
left = floor((length(Hlp_noise.Numerator)-1)/-2); %latency in samples
EDA_f = circshift(EDA_f, left); %correction
%replace last samples
EDA_f(length(EDA_f)+left:end) = EDA_f(length(EDA_f)+left-1);

%% Analyze for Artefacts

NaN_buffer_samp = floor(NaN_buffer.*SR/2); %NaN buffer in samples
art_buffer_samp = floor(art_buffer.*SR/2); %Artefact buffer in samples

EDA_f_NaN = EDA_f; %pre-allocate memory

% Difference between current sample and a art_buffer sample before
EDA_d = abs(EDA-circshift(EDA,-art_buffer_samp))*100; % x100 is for percent

EDA_d(end+left:end,1) = 0; %replace filter latency with zeros
EDA_d(1,1) = 0; %removes first change in case that EDA starts != 0

% Scan vector for artefacts and replace with NaNs

%while i < NaN_buffer_samp
for i = 1:NaN_buffer_samp
    if EDA_d(i,1)>thres
        EDA_f_NaN(1:i+NaN_buffer_samp) = NaN;
    end
end
%while NaN_buffer_samp<i<N-NaN_buffer_samp
for i = NaN_buffer_samp+1:N-NaN_buffer_samp
    if EDA_d(i,1)>thres
        EDA_f_NaN(i-NaN_buffer_samp:i+NaN_buffer_samp) = NaN;
    end
end
%while NaN_buffer_samp>N-NaN_buffer_samp
for i = N-NaN_buffer_samp+1:N
    if EDA_d(i,1)>thres
        EDA_f_NaN(i-NaN_buffer_samp:end) = NaN;
    end
end

% Scan vector for samples equal to range limits and replace with NaN
EDA_f_NaN(OutRange) = NaN;

%Evaluate quality of signal by counting NaNs (artefacts) in signal
Q = 100 - 100*(sum(isnan(EDA_f_NaN)))/length(EDA_f_NaN); %in percentage

%Evaluate difference before & after artefact, and add that factor to Q
%Search for first non NaN
i = 1;
while (isnan(EDA_f_NaN(i))==1 && i<N); i = i+1;end %Find 1st non NaN
%Calculate difference
while i<N
    if isnan(EDA_f_NaN(i))==0; i=i+1;
    else
        before = EDA_f_NaN(i-1); %record value before NaN
        j = i;
        while isnan(EDA_f_NaN(j))==1
            j=j+1;
            if j == N; break; end
        end
        after = EDA_f_NaN(j); %record value after NaN
        %Substract difference from Q (nansum is in case last sample==NaN)
        Q = Q - 50*abs(nansum([before -after])); % x50 is a manual scaling
        i = j+1;
    end
end

if Q <= 0; Q = 0; end;

if (debug==1);
fprintf('Confidence: %.2f%% - Filter length: %.2f[s]\n',Q,abs(left)/SR)
end

%% Interpolate between Artefacts

EDA_I = EDA_f_NaN;

%check if vector is only NaNs or no NaNs
if (sum(~(isnan(EDA_I)))>=2)&&(sum(isnan(EDA_I))>0)&&(InterOption==1)
    i = 1;
    while (isnan(EDA_I(i))==1 && i<N); i = i+1;end %Find 1st non NaN
    %replace 1st NaNs with 1st valid value
    if i~=1; EDA_I(1:i-1)=EDA_I(i); end
    
    while i<N
        if isnan(EDA_I(i))==0; i=i+1;
        else
            before = EDA_I(i-1); %record value before NaN
            j = i;
            while isnan(EDA_I(j))==1
                j=j+1;
                if j == N; EDA_I(i:N)=EDA_I(i-1); break; end
            end
            after = EDA_I(j); %record value after NaN
            EDA_I(i:j-1) = interp1([i j],[before after],i:j-1);
            i = j+1;
        end
    end
    
end

%Find new start value (could have changed due to artefacts)
EDA_o_T = EDA_o; %variable for Tonic offset
for Toffset=1:N
    if isnan(EDA_I(Toffset))==0 %find 1st non NaN
        EDA_o_T = EDA_o+EDA_I(Toffset);
        break;
    end
end

%% LP Filter (Tonic)

% Fc = [0.001 1];
% Astop = 60;        %stopband attenuation in dB
% Apass = 1;         %passband ripple in dB
% 
% % All frequency values are in Hz.
% Fpass = Fc(1);       % Passband Frequency
% Fstop = Fc(2);       % Stopband Frequency
% match = 'stopband';  % Band to match exactly
% 
% % Construct an FDESIGN object and call its BUTTER method.
% H2  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
% h2 = design(H2, 'butter', 'MatchExactly', match);

%Filter using loaded coefficients
EDA_T = filter(Hlp_tonic,EDA_I-EDA_I(Toffset));

%% HP Filter (phasic)

%Obtain phasic by substracting tonic from EDA_f
EDA_P = EDA_I-EDA_T;

%replace NaN using EDA_f_NaN (for no interpolation)
if InterOption~=1; EDA_P(isnan(EDA_f_NaN)) = NaN; end

%% Post-Processing

% Return vectors to original values in range and scale to match percent of
% range

EDA = (EDA+EDA_o)*100;
EDA_f = (EDA_f+EDA_o)*100;
EDA_f_NaN = (EDA_f_NaN+EDA_o)*100;
EDA_T = (EDA_T+EDA_o_T)*100;
EDA_P = (EDA_P+EDA_o)*100;
EDA_x = (EDA_I+EDA_o)*100; %Vector with NaNs or interpolation
EDA_I = (EDA_I+EDA_o)*100;

%% Plot

if debug == 1
    
    figure('name','EDA Tool','NumberTitle','off')
    
    %Create timeline
    t = (0:N-1).*1/SR;   
    x_limits = [t(1)-2 t(end)+2];
    
    %Y axis limits are min/max + 10% (unless range is < 30%)
    y_limits = [nanmin(nanmin([EDA_P EDA_T EDA_x]))*.9...
        nanmax(nanmax([EDA_P EDA_T EDA_x]))*1.1];
    if range(y_limits)<30
        y_limits=[mean(y_limits)-15 mean(y_limits)+15]; end
    if isnan(y_limits)==1; y_limits=0:1; end %In case all NaNs
    
    %vector to display filter latency
    latency = zeros(N,1)-1000; latency(N+left) = 2000;
    
    ax(1) = subplot(4,1,1);
    plot(t,EDA)
    title('Original EDA - Filtered EDA (EDA_f) - Filter Latency',...
        'Interpreter','none')
    hold on
    plot(t,EDA_f,'r','LineWidth',1)
    plot(t,latency,':','color','k'); %Shows filter latency
    hold off
    ylim([-0.05 100.05])
    xlim(x_limits)
    ylabel('range(%)')
    
    ax(2) = subplot(4,1,2);
    plot(t,EDA_d)
    hold on
    thres_plot = ones(N,1)*thres;
    plot(t,thres_plot,'k','LineWidth',2)
    ylim([0 thres*2.5])
    title(sprintf('Artefact Index - Window %.2f[s], Threshold %.1f%%',...
        art_buffer,thres))
    xlim(x_limits)
    ylabel('range(%)')
    
    ax(3) = subplot(4,1,3);
    plot(t,EDA_f_NaN,'k','LineWidth',1)
    hold on
    plot(t(isnan(EDA_f_NaN)),EDA_I(isnan(EDA_f_NaN)),'r',...
    'LineStyle', 'none', 'Marker', '.', 'MarkerSize',6)
    ylabel('range(%)')
    ylim(y_limits)
    title(sprintf('EDA without artifacts (EDA_x) - Confidence: %.1f%%',...
        Q),'Interpreter','none')
    xlim(x_limits)
    
    ax(4) = subplot(4,1,4);
    plot(t,EDA_T,'r','LineWidth',1)
    title('Phasic EDA (EDA_P) - Tonic EDA (EDA_T)','Interpreter','none')
    hold on
    plot(t,EDA_P,'b','LineWidth',1)
    hold off
    ylim(y_limits)
    xlim(x_limits)
    xlabel('time(s)')
    ylabel('range(%)')
    linkaxes(ax,'x')
    legend('Tonic','Phasic','Location','NorthWest')
end

%% Resample

if ResampleOption == 1    
    t_o = (0:N_o-1)./SR_original; %time vector for original signal    
    EDA_f = interp1(t,EDA_f,t_o); %resample EDA to original SR
    EDA_x = interp1(t,EDA_x,t_o); %this may cause warnings with NaNs
    EDA_P = interp1(t,EDA_P,t_o);
    EDA_T = interp1(t,EDA_T,t_o);
end

end