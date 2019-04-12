%----------------------------------------------------------------------
%BrainAmpConverter v1.7
%
%Conversion of physiological data recorded during fMRI scanning using
%Brainamp systems to HeRa (heart rate) and Autonomate (SCR) formatted
%files.
%
%Requires Fieldtrip IO functions for reading the EEG file created by Brainamp.
%These are included in the fileio subfolder
%
%v1.0 September 2015; changes wrt earlier scripts incorporated into this:
%- All previous conversion scripts were combined into brainampconverter
%- Added detection of fMRI scan runs and separate conversion per run
%- Added a deconvolution filter to remove scanner artifact
%- Switched to wider canonicals to take context into account
%- Changed peak selection algorithm from gaussian mixture
%model to a method that minizes the imbalance between positive extremes
%and negative extremes in the IBI timeseries.
%
%v1.1 (20151006 for hanjian)
%- added a third channel for Autonomate (with all event codes set to 1)
%- moved channel names to the settings
%
%v1.2 (20151009 for Marloes/Piet's analyses of BIG-C data
%- Improved channel name recognition (multiple options in settings)
%- Fix for data acquisition that to briefly before scan starts.
%- Now converts the entire file in addition to every fMRI sequence found.
%- Added Stimulus events to the HeRa file.
%- Added respiration data to the SCR file for autonomate (for future use)
%- Moved defaults to separate m-file
%- Fixed bug in canonical pulse window setting
%
%v1.3 Fixed bug in temporal delay in pulse data after filtering (20151026).
%
%v1.4 Fixed bug in discarding the first trigger pulse of the first fMRI
%     run that is detected. (20151119)
%v1.5 Updated the filio folder from FieldTrip to avoid crashes on newer
%     Matlab versions. (20180530)
%v1.6 Fixed problem with downsampling of SCR data when the downsample
%     factor is not an integer.
%v1.7 Fixed problem in some files in which the outlier minimization
%     procedure resulted in a wrong threshold, removing all peaks.
%EJH 2011-18
%----------------------------------------------------------------------
function brainampconverter(eegFilename)

if ~nargin
    [eegfile, eegpath, irr] = uigetfile( ...
       {'*.eeg','EEG-files (*.eeg)'}, ...
        'Pick an EEG file');
    eegFilename = fullfile(eegpath,eegfile);
    if numel(eegpath)==1; return; end
    
end

%Get default settings
brainampconverter_defaults;                             

%Add fieldtrip's IO functions to the path
[cpath,irr]=fileparts(mfilename('fullpath'));
addpath(fullfile(cpath,'fileio'));

%Display loading
disp(['... Loading ',eegFilename])

%Read the brainamp files
hdr=ft_read_header(eegFilename);
dat=ft_read_data(eegFilename);
evt=ft_read_event(eegFilename);

%Calculate downsample factors
HRdsFactor = hdr.Fs/HRdsFreq;
SCRdsFactor = hdr.Fs/SCRdsFreq;

%Search for the SCR, HR, and Resp channels
HRchannel = [];
Respchannel = [];
SCRchannel = [];
chanRec = zeros(1,numel(hdr.label));
for i=1:numel(hdr.label)
    
    %If it's an SCR channel
    if ismember(hdr.label{i},SCRchannelnames)
        SCRchannel=i;
        chanRec(i)=1;
        disp('... Found skin conductance recording channel')
    elseif ismember(hdr.label{i},HRchannelnames) %If HR
        HRchannel=i;
        chanRec(i)=1;
        disp('... Found finger pulse HR recording channel')
    elseif ismember(hdr.label{i},Respchannelnames) %if Respiration
        Respchannel=i;
        chanRec(i)=1;
        disp('... Found respiration recording channel')
    end
end

%Find unrecognized channels, if any
unrecChan = find(chanRec==0);

%Give warning if there were unrecognized channels
for i=1:numel(unrecChan)
    disp(['Warning: Unable to identify brainamp channel named "',hdr.label{i},'" [in hdr.label{',num2str(i),'}].']);
end

%Extract all scan triggers
cell_evt=permute(struct2cell(evt),[3,1,2]); %Restructure the data for easier searching

%Try to figure out what code was used for triggers
for i=1:numel(ResponseCodes)
    rCinds{i} = find(strcmp(cell_evt(:,1),ResponseCodes{i})); %Get indices of response codes
    nrCinds(i) = numel(rCinds{i});
end
%Figure out which one is correct; if one is at least three,
%and all others are less than three
if sum(nrCinds>2)==1 %only one bigger than 2
    %Get all the response codes
    [irr,cRCi]=max(nrCinds);
    ResponseCode = ResponseCodes{cRCi};
    rCinds = find(strcmp(cell_evt(:,1),ResponseCode)); %Get indices of response codes
    allResponseCodes=cell2mat(cell_evt(rCinds,2)); %Get all trigger/response codes
    allResponseCodes(double(allResponseCodes)<48 |double(allResponseCodes)>57)=' '; %remove text from strings
    allResponseCodesNum=str2num(allResponseCodes); %Convert to numbers
    triggerLines = rCinds(bitget(allResponseCodesNum,TriggerBit)==1); %lines in cell_evt with triggers
    triggerSamps = cell2mat(cell_evt(triggerLines,3)); %Get time stamps of all trigger codes
    triggerSamps=triggerSamps([1;find(diff(triggerSamps)>(hdr.Fs/10))+1]); %Remove double triggers (occurring within 100ms, due to overlapping response codes)

else
    %Continue without response codes (therefore no triggers and no fMRI runs)
    disp('... Could not find any Response codes (needed to detect scanner triggers). Please check the "ResponseCodes" setting in brainampconverter. Converting only full file.');
    ResponseCode = [];
    allResponseCodes = [];
    allResponseCodesNum = [];
    triggerLines = [];
    triggerSamps = [];
end

%Figure out what is the correct stimulus code
for i=1:numel(StimulusCodes)
    sCinds{i} = find(strcmp(cell_evt(:,1),StimulusCodes{i})); %Get indices of stimulus codes
    nsCinds(i) = numel(sCinds{i});
end

%Figure out which one is correct; there should be only one present,
%otherwise display a warning and continue without
if sum(nsCinds>0)==1 %only one bigger than 2
    [irr,cSCi]=max(nsCinds);
    StimulusCode = StimulusCodes{cSCi};
    %Extract all stimulus codes
    sCinds = find(strcmp(cell_evt(:,1),StimulusCode)); %Get indices of stimulus codes
    allStimulusCodes=cell2mat(cell_evt(sCinds,2)); %Get all stimulus codes
    allStimulusCodes(double(allStimulusCodes)<48 |double(allStimulusCodes)>57)=' '; %remove text from strings
    allStimulusCodesNum=str2num(allStimulusCodes); %Convert to numbers
    stimulusSamps = cell2mat(cell_evt(sCinds,3)); %Get time stamps of all trigger codes
    disp(['... Found ',num2str(numel(allStimulusCodesNum)),' stimulus code events.'])
elseif sum(nsCinds>0)>1
   disp('Cannot determine which code is used for stimuli. Please check the StimulusCodes setting in brainampconverter. Continuing without stimulus event markers.'); 
   allStimulusCodesNum = [];
   stimulusSamps = [];
elseif sum(nsCinds>0)==0
   disp('Cannot find stimulus codes. Please check the StimulusCodes setting in brainampconverter. Continuing without stimulus event markers.'); 
   allStimulusCodesNum = [];
   stimulusSamps = [];
end




%If scan triggers have been found, try to detect regularities to find
%MRI scan runs
runTriggers=[];
runStimulusSamps=[];
runStimulusCodes=[];
if numel(triggerLines)>0
    %Detect starts and ends of runs (at least two TRs)
    samp10ms = hdr.Fs/100; %How many samples is 10 ms?
    runOnsets = find(diff(vertcat(0,abs(diff(diff(triggerSamps)))<samp10ms,0))==1);
    runOffsets = find(diff(vertcat(0,abs(diff(diff(triggerSamps)))<samp10ms,0))==-1)+1;
    %Figure out if run has been detected
    if numel(runOnsets)==0 | numel(runOnsets)~=numel(runOffsets)
       disp('... Cannot figure out fMRI runs. Converting only full file.');
    else %if at least one fMRI run was detected
        runTriggers = [];
        for i=1:numel(runOnsets)
            %Get the time stamps of triggers for each run
            runTriggers{i} = triggerSamps(runOnsets(i):runOffsets(i));
            %Get the time stamps of stimuli for each run
            runStimulusSamps{i} = stimulusSamps(stimulusSamps>=runTriggers{i}(1) & stimulusSamps<=runTriggers{i}(end));
            %Get the stimulus codes for each run
            runStimulusCodes{i} = allStimulusCodesNum(stimulusSamps>=runTriggers{i}(1) & stimulusSamps<=runTriggers{i}(end));
        end
        disp(['... Found ',num2str(numel(runTriggers)),' fMRI run(s)'])
    end
end


%Now loop over conversions, starting with full file, and then
%separately for fMRI runs
for i=(1-FullFileConversion):numel(runTriggers)

    %Is this a MRI run conversion?
    if i==0 %convert entire file
        cMRIconv = false;
        cextraWindow = 0; %There is no extra data for the entire file
    else %convert fMRI run
        cMRIconv = true;
        %Calculate the approximate TR in ms for this run
        aTR=(mean(diff(runTriggers{i}))/hdr.Fs)*1000;
        cextraWindow = extraWindow; %Cut out eg 10s more before and after (see settings)
    end
    
 
    %-----------------------------------------------------------------
    %If there is heart rate data, convert it to Hera
    %If present, respiration data is included in this file
    %-----------------------------------------------------------------
    
    if numel(HRchannel)>0
        
        %Pre-specify the low-pass filter for HR (because time shift needs to be
        %known in advance)
        Fst = 3; %Fst ? frequency at the end of the stop band. was 3 - 4
        Fp = 4; % Fp ? frequency at the start of the pass band.
        Ap =.5; %Ap ? amount of ripple allowed in the pass band in decibels
        Ast=40; %Ast ? attenuation in the stop band in decibels (the default units)
        d_lpf = fdesign.lowpass('Fp,Fst,Ap,Ast',Fst,Fp,Ap,Ast,(hdr.Fs/HRdsFactor));
        Hd_lpf = design(d_lpf,'equiripple');
        phaseShiftSamples_lpf = numel(Hd_lpf.Numerator)/2;

        %Pre-specify the high-pass filter for HR
        cFst = .2; %Fst ? frequency at the end of the stop band. Specified in normalized frequency units. Also called Fstop.
        cFp = .3; % Fp ? frequency at the start of the pass band. Specified in normalized frequency units. Also called Fpass.
        rcFst = (2*pi)/((1/cFst)*(hdr.Fs/HRdsFactor));
        rcFp = (2*pi)/((1/cFp)*(hdr.Fs/HRdsFactor));
        Ast = 60; %Ast ? attenuation in the stop band in decibels (the default units). Also called Astop.
        Ap = 1;   %Ap ? amount of ripple allowed in the pass band in decibels (the default units). Also called Apass.
        d_hpf = fdesign.highpass('Fst,Fp,Ast,Ap',rcFst,rcFp,Ast,Ap);
        Hd_hpf = design(d_hpf,'equiripple');
        phaseShiftSamples_hpf = numel(Hd_hpf.Numerator)/2;
        %Calculate total filter shift
        phaseShiftSamples=round(phaseShiftSamples_lpf+phaseShiftSamples_hpf);

        %Get start of run - extra window - or set to 1 for entire file
        if cMRIconv
            startRunWin = round(runTriggers{i}(1) - (cextraWindow*hdr.Fs));

            %Get end of run + 1 TR + extra window + filter shift for HR
            endRunWin = round(runTriggers{i}(end) + ...
                (cextraWindow*hdr.Fs) + (aTR/1000)*hdr.Fs + ...
                phaseShiftSamples*HRdsFactor);
            
            %Correct timing of scan triggers
            triggerSamps = runTriggers{i} - (startRunWin-1);
            
            %If there are stimulus events in this run, correct their timing
            if numel(runStimulusSamps)>0 && numel(runStimulusSamps{i}>0)
                crunStimulusSamps = runStimulusSamps{i} - (startRunWin-1); %do these exist?
                crunStimulusCodes = runStimulusCodes{i};
            else
                crunStimulusSamps = [];
                crunStimulusCodes = [];
            end
            
        else %Set to entire file 
            startRunWin = 1;
            endRunWin = size(dat,2);
            crunStimulusSamps = stimulusSamps;          %Copy all entire stimulus events
            crunStimulusCodes = allStimulusCodesNum;    %and their codes

        end
        endRunWin = min([endRunWin,size(dat,2)]); %cut off if too long

        %Extract the pulse/resp data
        respdat=[];
        if startRunWin>0 %If the extrawindow does not set the start point before time point 1
            pulsedat = dat(HRchannel,startRunWin:endRunWin)';
            if numel(Respchannel)>0
                respdat = dat(Respchannel,startRunWin:endRunWin)'; %should not crash even if no resp data
            end
        else
            pulsedat = [zeros(1,(-startRunWin+1)), ... %Put zeros before
                dat(HRchannel,1:endRunWin)]';
            if numel(Respchannel)>0
                respdat = [zeros(1,(-startRunWin+1)), ... %same for resp
                    dat(Respchannel,1:endRunWin)]';
            end
        end
        
        %Downsample the data
        dspulsedat = downsample(pulsedat,HRdsFactor);
        dsrespdat = downsample(respdat,HRdsFactor);  
        
        %Keep the raw data (for debugging)
        %dspulsedat_raw = dspulsedat;

        %Downsample the triggers
        dstriggerSamps = round(triggerSamps/HRdsFactor);
        dscrunStimulusSamps = round(crunStimulusSamps/HRdsFactor);
        
        %If this is an fMRI run, remove scanner artifacts
        if cMRIconv
            %Estimate the shape of the TR-onset-timelocked signal modulation
            disp(['... Filtering out scanner artifacts in HR data for run ',num2str(i)])
            TRlocked=[];
            for j=1:numel(dstriggerSamps)
                TRlocked(j,:) = dspulsedat(round(dstriggerSamps(j)-((diff(dstriggerSamps))*.5)):...
                    round(dstriggerSamps(j)+((diff(dstriggerSamps))*.5)))';
            end

            %Standardize this
            meanTRlocked=mean(TRlocked);
            meanTRlocked=meanTRlocked-mean(meanTRlocked);
            meanTRlocked=meanTRlocked./std(meanTRlocked);

            %Create deconvolution model
            dm=zeros(numel(dspulsedat),numel(dstriggerSamps));
            for j=1:numel(dstriggerSamps)
                dm(round(dstriggerSamps(j)-((diff(dstriggerSamps))*.5)):...
                    round(dstriggerSamps(j)+((diff(dstriggerSamps))*.5)),j) = ...
                    meanTRlocked';
            end

            %Deconvolve the data / keep residual
            b = inv(dm'*dm)*dm'*dspulsedat;
            dspulsedat = dspulsedat - dm*b; 
        end %scanner artifact correction
        
        %Apply the low-pass filter
        lpf_dspulsedat = filter(Hd_lpf,dspulsedat);

        %Apply the high-pass filter
        hpf_lpf_dspulsedat = filter(Hd_hpf,lpf_dspulsedat);  

        %Shift the signal back to correct for filter-induced time delay
        ts_hpf_lpf_dspulsedat=...
            hpf_lpf_dspulsedat(1+phaseShiftSamples:end);
        %Then cut respiration to the same length (if present)
        if numel(dsrespdat)>0
            dsrespdat=dsrespdat(1:numel(ts_hpf_lpf_dspulsedat));
        end
            
        %Check if converted puls file already exists
        [pulspath,pulsfile]=fileparts(eegFilename);

        %Create the ID for the filename
        if cMRIconv
            runID = ['_fMRI_run_',num2str(i),'_',...
                'TR',num2str(round(aTR)),'ms_',...
                num2str(numel(triggerSamps)),'vols_',...
                num2str(round(startRunWin/hdr.Fs)),'-',...
                num2str(round((endRunWin/hdr.Fs)-(phaseShiftSamples/(hdr.Fs/HRdsFactor)))),'s'];
        else
            runID = ['_full_recording_',...
                num2str(round(numel(ts_hpf_lpf_dspulsedat)/(hdr.Fs/HRdsFactor))),'s'];
        end
        
        pulsfilepath=fullfile(pulspath,[pulsfile,runID,'_hera.puls']);
        pulsmatfilepath=fullfile(pulspath,[pulsfile,runID,'_hera.mat']);

        %If this file exists already, do not overwrite because
        %changes to the file might be lost
        if exist(pulsfilepath)&avoidOverwrite

            continueoverwr = questdlg(['Warning: Converted Hera Puls file ',pulsfilepath,...
                ' already exists and may have been edited. Continue to overwrite this file?'],...
                'Overwrite file?', 'No');
            if  strcmp(continueoverwr,'No') | strcmp(continueoverwr,'Cancel')
                error('Stopped conversion.')
            end
        end

        %Run the conversion
        %Run the PPU pulse detection
        if cMRIconv
            disp(['... Running pulse peak detection and conversion to HeRa for ',eegFilename,' fMRI run #',num2str(i)])
        else
            disp('... Running pulse peak detection and conversion to HeRa for full recording')
        end
        
        %Note: adjusted use of time shifted pulse (20151026)
        [ccspeaks,ibits]=brainampconverter_pulsepeak(ts_hpf_lpf_dspulsedat,hdr.Fs/HRdsFactor,canpulsewin,visualize);

        %Make a puls/mat file pair readable by HeRa
        matfile.settings= struct(...
            'pulsefile_extension', 'puls', ...
            'prepeak', 5000, ...
            'prepeakdelay', 2,...
            'samplerate', hdr.Fs/HRdsFactor, ...
            'barscalefactor', 0.57, ...
            'hpd', -1, ...
            'lpd', -1, ...
            'avgHeartbeat', 1, ...
            'ballmobilityrange', 0.1, ...
            'startCutout', 2, ...
            'endCutout', 2, ...
            'hpt', 0.0050, ...
            'lpt', 1, ...
            'lowerboundLF', 0.0500, ...
            'upperboundLF', 0.1500, ...
            'lowerboundHF', 0.1500, ...
            'upperboundHF', 0.4000);

        %Set the raw pulse data to mean 2000 and SD 500
        %To make it correspond to Siemens pulse data for HeRa

        %Note: corrected time shift pulse (20151026)
        matfile.bpf_pulsedata = ts_hpf_lpf_dspulsedat'; %Needs to be horizontal for Hera
        %Cut off raw pulse data to the same length if longer (20151026)
        matfile.rawpulsedata = (round(2000+((dspulsedat-mean(dspulsedat))./std(dspulsedat)).*500))'; %Needs to be horizontal for Hera
        matfile.rawpulsedata = matfile.rawpulsedata(1:numel(matfile.bpf_pulsedata));
        
        matfile.prepeaklocs = ccspeaks';
        matfile.prepeakTimes = (matfile.prepeaklocs-1).*(HRdsFactor/hdr.Fs);
        matfile.preibitimeseries = diff(matfile.prepeakTimes);
        matfile.preibichannel = (ibits.*(1000/(hdr.Fs/HRdsFactor)))';
        matfile.prereject = {};
        matfile.prepulsezoom = [0 1];
        matfile.postibichannel = [];
        matfile.sampleTime = (0:numel(pulsedat)-1)./(hdr.Fs/HRdsFactor);
        matfile.scalFreq=[];
        matfile.aBin=[];
        matfile.CWT=[];
        matfile.topTrace=[];
        matfile.spectrumFrequencies=[];
        matfile.spectrumPower=[];
        matfile.lowerboundLFIndex=[];
        matfile.upperboundLFIndex=[];
        matfile.lowerboundHFIndex=[];
        matfile.upperboundHFIndex=[];
        matfile.outmeasures.preBPM = -1;
        matfile.outmeasures.rMSSD = -1;
        matfile.outmeasures.LFPower =[];
        matfile.outmeasures.HFPower =[];
        matfile.outmeasures.LFHFratio =[];
        matfile.outmeasures.postBPM =[];

        %Add triggers as markers in the hera file
        matfile.markerlocs = dstriggerSamps;
        
        %Add stimulus codes
        matfile.StimulusSamps = dscrunStimulusSamps;
        matfile.StimulusCodes = crunStimulusCodes;

        %Include respiration data (although not used by HeRa)
        matfile.rawrespdata = dsrespdat';

        %Save the files; first an empty (used only for identification) puls file
        fid=fopen(pulsfilepath,'w');
        fclose(fid);

        save(pulsmatfilepath,'matfile');
    end
    
    
    %-----------------------------------------------------------------
    %If there is SCR data, convert it to ASCII for use with Autonomate
    %-----------------------------------------------------------------

    if numel(SCRchannel)>0
        
        %Check what is the mean SCR level, and warn if out of range.
        if ~cMRIconv
            meanSCR = mean(dat(SCRchannel,:))/uVperuSiem;
            if meanSCR<1
                disp(['Warning: Please note that mean skin conductance level is low: ',...
                    num2str(meanSCR),' microSiemens. Make sure the calibration setting (uVperuSiem) in brainampconverter_defaults.m is correct.'])
            elseif meanSCR>30
                disp(['Warning: Please note that mean skin conductance level is high: ',...
                    num2str(meanSCR),' microSiemens. Make sure the calibration setting (uVperuSiem) in brainampconverter_defaults.m is correct.'])
            end
        end
        
        %Pre-specify the low-pass filter for SCR (because time shift needs to be
        %known in advance)
        Fst = 1.5; %Fst ? frequency at the end of the stop band. was 3 - 4
        Fp = 2; % Fp ? frequency at the start of the pass band.
        Ap =.5; %Ap ? amount of ripple allowed in the pass band in decibels
        Ast=40; %Ast ? attenuation in the stop band in decibels (the default units)
        d_lpf = fdesign.lowpass('Fp,Fst,Ap,Ast',Fst,Fp,Ap,Ast,(hdr.Fs/HRdsFactor));
        Hd_lpf = design(d_lpf,'equiripple');
        phaseShiftSamples_lpf = round(numel(Hd_lpf.Numerator)/2);
        
        %Pre-specify moving average filter for SCR (for time shift)
        %LPF + moving average is in accordance with recommendation for
        %autonomate (Green et al 2013)
        nSamplesmaf = SCR_MAFwidth*SCRdsFreq;
        maf = (1/nSamplesmaf)*ones(nSamplesmaf,1);
        phaseShiftSamples_maf = nSamplesmaf/2;
        
        %Calculate total time shift
        phaseShiftSamples = round(phaseShiftSamples_lpf+phaseShiftSamples_maf);
        
        %Get start of run - extra window for SCR - or set to 1 for entire file
        if cMRIconv
            startRunWin = round(runTriggers{i}(1) - (cextraWindow*hdr.Fs));
            
            %Get end of run + 1 TR + extra window + filter shift for SCR
            endRunWin = round(runTriggers{i}(end) + ...
                (cextraWindow*hdr.Fs) + (aTR/1000)*hdr.Fs + ...
                phaseShiftSamples*SCRdsFactor);
            endRunWin = min([endRunWin,size(dat,2)]); %cut off if too long

            %Correct timing of scan triggers and stimulus codes
            triggerSamps = runTriggers{i} - (startRunWin-1);
            
            %If there are stimulus events in this run, correct their timing
            if numel(runStimulusSamps)>0 && numel(runStimulusSamps{i}>0)
                crunStimulusSamps = runStimulusSamps{i} - (startRunWin-1); %do these exist?
                crunStimulusCodes = runStimulusCodes{i};
            else
                crunStimulusSamps = [];
                crunStimulusCodes = [];
            end
            
        else %Set to entire file 
            
            startRunWin = 1;
            endRunWin = size(dat,2);
            crunStimulusSamps = stimulusSamps;          %Copy all entire stimulus events
            crunStimulusCodes = allStimulusCodesNum;    %and their codes

        end
        endRunWin = min([endRunWin,size(dat,2)]); %cut off if too long
        
        %Extract the SCR data, include respiration if possible
        respdat=[];
        if startRunWin>0 %If the extrawindow does not set the start point before time point 1
            SCRdat=dat(SCRchannel,startRunWin:endRunWin)';
            if numel(Respchannel)>0
                respdat = dat(Respchannel,startRunWin:endRunWin)'; %should not crash even if no resp data
            end        
        else
            SCRdat = [zeros(1,(-startRunWin+1)), ... %Put zeros before
                dat(SCRchannel,1:endRunWin)]';
            if numel(Respchannel)>0
                respdat = [zeros(1,(-startRunWin+1)), ... %same for resp
                    dat(Respchannel,1:endRunWin)]';
            end
        end
        
        %Downsample the data
        %20180530: when the required downsample factor is not an integer,
        %the downsample function crashes. For backward compatibility, the
        %old downsampling is still used if the downsample factor is an
        %integer, otherwise the data is resampled at the new frequency
        
        if floor(SCRdsFactor)==SCRdsFactor %if integer
            dsSCRdat = downsample(SCRdat,SCRdsFactor);
            dsrespdat = downsample(respdat,SCRdsFactor); 
        else %if not integer
            dsSCRdat = resample(SCRdat,SCRdsFreq,hdr.Fs);
            if length(respdat)>0 %if there's respiration data, do the same
                dsrespdat = resample(respdat,SCRdsFreq,hdr.Fs); 
            else
                dsrespdat = [];
            end
        end            
       
        %Keep the raw data (for debugging)
        %dsSCRdat_raw = dsSCRdat;

        %Downsample the triggers / stimulus events
        dstriggerSamps = round(triggerSamps/SCRdsFactor);
        dscrunStimulusSamps = round(crunStimulusSamps/SCRdsFactor);

        %If this is an fMRI run, remove scanner artifacts
        if cMRIconv
        
            %Estimate the shape of the TR-onset-timelocked signal modulation
            disp(['... Filtering out scanner artifacts in SCR data for run ',num2str(i)])
            TRlocked=[];
            for j=1:numel(dstriggerSamps)
                TRlocked(j,:) = dsSCRdat(round(dstriggerSamps(j)-((diff(dstriggerSamps))*.5)):...
                    round(dstriggerSamps(j)+((diff(dstriggerSamps))*.5)))';
            end

            %Standardize this
            meanTRlocked=mean(TRlocked);
            meanTRlocked=meanTRlocked-mean(meanTRlocked);
            meanTRlocked=meanTRlocked./std(meanTRlocked);

            %Create deconvolution model
            dm=zeros(numel(dsSCRdat),numel(dstriggerSamps));
            for j=1:numel(dstriggerSamps)
                dm(round(dstriggerSamps(j)-((diff(dstriggerSamps))*.5)):...
                    round(dstriggerSamps(j)+((diff(dstriggerSamps))*.5)),j) = ...
                    meanTRlocked';
            end

            %Deconcolve the data / keep residual
            b = inv(dm'*dm)*dm'*dsSCRdat;
            ac_dsSCRdat = dsSCRdat - dm*b; 
        else
            ac_dsSCRdat = dsSCRdat;
        end
        
        %Apply low pass filter
        lpf_ac_dsSCRdat = filter(Hd_lpf,ac_dsSCRdat);

        %Apply moving averate filter
        maf_lpf_ac_dsSCRdat = filter(maf,1,lpf_ac_dsSCRdat);

        %Shift the signal back to correct for filter-induced time delay
        ts_maf_lpf_ac_dsSCRdat = ...
            maf_lpf_ac_dsSCRdat(1+phaseShiftSamples:end);
        %Then cut respiration to the same length (if present)
        if numel(dsrespdat)>0
            dsrespdat=dsrespdat(1:numel(ts_maf_lpf_ac_dsSCRdat));
        end

        %Create the Autonomate-compatible event channel
        Aut_EventChannel = zeros(numel(ts_maf_lpf_ac_dsSCRdat),1);
        Aut_EventChannel(dscrunStimulusSamps)=crunStimulusCodes;
        
        %And add another column containing only one event type
        Aut_EventChannel_bool = Aut_EventChannel~=0;
        
        %Scale exported SCR data to microSiemens (see header above)
        ts_maf_lpf_ac_dsSCRdat = ts_maf_lpf_ac_dsSCRdat./uVperuSiem;
        
        %Export the data
        %Create the ID for the filename
        if cMRIconv
            runID = ['_fMRI_run_',num2str(i),'_',...
                'TR',num2str(round(aTR)),'ms_',...
                num2str(numel(triggerSamps)),'vols_',...
                num2str(round(startRunWin/hdr.Fs)),'-',...
                num2str(round((endRunWin/hdr.Fs)-(phaseShiftSamples/(hdr.Fs/SCRdsFactor)))),'s'];
                disp(['... Running conversion of SCR to Autonomate for ',eegFilename,' fMRI run #',num2str(i)])

        else
            runID = ['_full_recording_',...
                num2str(round(numel(ts_maf_lpf_ac_dsSCRdat)/(hdr.Fs/SCRdsFactor))),'s'];
                disp('... Running conversion of SCR to Autonomate for full recording')
        end
        
        [autpath,autfile]=fileparts(eegFilename);
        autfilepath=fullfile(autpath,[autfile,runID,'_autonomate.txt']);
        if numel(dsrespdat)>0
            outSCR = [ts_maf_lpf_ac_dsSCRdat,dsrespdat,Aut_EventChannel,Aut_EventChannel_bool];
        else
            outSCR = [ts_maf_lpf_ac_dsSCRdat,Aut_EventChannel,Aut_EventChannel_bool];
        end
        
        dlmwrite(autfilepath,outSCR,' ');
    end
end




%--------------------------------------------------------------------------
%brainampconverter_pulsepeak
%
%Detects peaks in pulse oximeter data from BIOPAC or BRAINAMP systems
%
%Using these steps:
%1. Standardize the pulse waveform
%2. Determine all peaks present in the signal (candidate heart beats)
%3. Get all certain peaks with a very stringent threshold
%4. Determine a "canonical" pulse waveform from subset of certain peaks
%5. Calculate correlation (/ fisher z) of local signal around each peak in the signal
%   with this canonical waveform.
%6. Estimate the overall BPM by calculating the median BPM as a function of z threshold.
%7. Loop again over all thresholds to find at which threshold 
%   you get approximately equally many upward as downward extremes (relative to the estimaes BPM)
%   in the IBI timeseries
%8. Detect and try to repair IBIs that are too short
%9. Detect and try to repair IBIs that are too long
%10.Calculate new IBI timeseries and visualize results
%
%EJH 2011-15
%--------------------------------------------------------------------------
function [speaks,ibits]=brainampconverter_pulsepeak(fpulsedat,sR,canpulsewin,visualize)

%Settings
nSD = 2;              %Set number of SDs deviation for initial "certain" peak detection
minapeakdistms = 250; %minimum distance between two detected peaks
outlierIBIth = 1.75;  %Threshold when counting outlier IBIs to determine threshold
certainPeakRange=[11:30]; %if 11;30: take the 11th strongest peak until 30th as certain peaks because strongest are often artefact
IBIsdco = 2; %3         %Setting for step 8/9 (remove too short/long IBIs): IBI cutoff for IBIs that are too short/long
nSDrep = 1.5;         %Setting for step 9 (remove too long IBIs): Criterion for standard deviation of repair peak
nSDtrue = 4;          %Setting for step 9 (remove too long IBIs): Maximum number of SD from mean of accepted peaks distribution for repair peaks

%--------------------------------------------------------------------------
%Step 1: standardize BPM channel
sfpulsedat=(fpulsedat-mean(fpulsedat))./std(fpulsedat);

%--------------------------------------------------------------------------
%Step 2: Determine all peaks
apeaks = find(diff(sign(diff(sfpulsedat)))<0)+1;

%If two adjacent peaks are closer than minapeakdistms (eg 250 ms), then take only the
%highest of the two
while sum(diff(apeaks)<(minapeakdistms/1000)*sR) > 0
    %Find the smallest gap between peaks
    [irr,ind_min]=min(diff(apeaks));
    selvec = true(size(apeaks));
    selvec(ind_min)=false;
    apeaks=apeaks(selvec);
end

%Maria's method instead:
%apeaks = peakfinder(sfpulsedat);
%--------------------------------------------------------------------------
%Step 3: Determine only absolutely certain peaks
%Take the range defined in certainPeakRange above
%Because among the strongest 10 are often artifacts

strongPeaks = flipud(sortrows([sfpulsedat(apeaks),[1:numel(apeaks)]']));
cpeaks = apeaks(strongPeaks(certainPeakRange,2));
cpeaks = sortrows(cpeaks); %order them chrnonologically

%--------------------------------------------------------------------------
%Step 4 determine the canonical pulse (or pulses)
cphwin=round((canpulsewin*sR)/2);
canpulse=zeros(cphwin*2,1);
ilist=find(cpeaks>(.5*canpulsewin*sR) & cpeaks<(numel(sfpulsedat)-.5*canpulsewin*sR))';
for i=ilist
    canpulse=canpulse+sfpulsedat(cpeaks(i)-cphwin+1:cpeaks(i)+cphwin);
end
canpulse=canpulse./numel(ilist);

if visualize
    figure
    plot(canpulse)
    title 'Canonical Pulse'
end

%--------------------------------------------------------------------------
%Step 5: Loop over all peaks to determine their shape correspondence with
%canonical pulse beat and calculate fisher transformed z
fzapeaks=nan(numel(apeaks),1);
for i=find(apeaks>cphwin &apeaks<numel(sfpulsedat)-cphwin)'
    cc=corrcoef(sfpulsedat(apeaks(i)-cphwin+1:apeaks(i)+cphwin),canpulse);
    fzapeaks(i)=.5*(log(1+cc(2,1)) - log(1-cc(2,1)));
end

if visualize
    figure
    cla
    hold on
    plot(sfpulsedat)
    vis_fzapeaks = zeros(size(sfpulsedat));
    vis_fzapeaks(apeaks) = fzapeaks;
    vis_fzapeaks=(vis_fzapeaks.*10)-20;
    plot(vis_fzapeaks,'k')
end



%--------------------------------------------------------------------------
%Step 6: estimate the BPM by calculating the median BPM as a function of
%threshold, and then taking the mean of those

j=1;
xax=min(fzapeaks):.01:max(fzapeaks);
for i=xax
    th_apeaks = apeaks(fzapeaks>i);
    th_ibits = (diff(th_apeaks)./sR)*1000;
    th_medbpm(j) = 60000/median(th_ibits);
    j=j+1;
end

%Estimate the median heart rate frequency
medbpm = nanmedian(th_medbpm);

%--------------------------------------------------------------------------
%Step 7: Loop again over all thresholds to find at which threshold 
%you get approximately equally many upward as downward extremes (relative to the estimaes BPM)
%in the IBI timeseries

j=1;
xax=min(fzapeaks):.01:max(fzapeaks);

for i=xax
    th_apeaks = apeaks(fzapeaks>i);
    th_ibits = (diff(th_apeaks)./sR)*1000;
    %Calculate deviations
    meddevs = th_ibits - (60000/medbpm);
    %calculate outlier thresholds
    medIBI = 60000/medbpm;
    cmaxIBI = (medIBI*outlierIBIth)-medIBI;
    cminIBI = (medIBI/outlierIBIth)-medIBI;
    %Count when exceeding IBI outlier threshold
    posmeddevs(j) = sum(meddevs> cmaxIBI) / numel(th_ibits);
    negmeddevs(j) = sum(meddevs< cminIBI) / numel(th_ibits);
    j=j+1;
end

%Calculate the threshold applied to the z scores by finding the threshold
%at which there are equally many positive as negative deflections

%20181019: force the chosen value to be be between the two peaks (of
%posmeddevs and negmeddevs)
[a,peakloc_negmeddevs]=max(negmeddevs);
peakloc_negmeddevs=peakloc_negmeddevs(end);
[a,peakloc_posmeddevs]=max(posmeddevs);
peakloc_posmeddevs=peakloc_posmeddevs(1);
%Error if peakloc_negmeddevs is not before peakloc_posmeddevs
if peakloc_posmeddevs<=peakloc_negmeddevs
    error('Cannot find optimal cut-off for minimizing outliers in peak detection');
end

%cut this part out for both
cut_negmeddevs=negmeddevs(peakloc_negmeddevs:peakloc_posmeddevs);
cut_posmeddevs=posmeddevs(peakloc_negmeddevs:peakloc_posmeddevs);

%Find the point where these two come closest
[a,b]=min(abs(cut_posmeddevs-cut_negmeddevs));
%Correct b to account for removed part before peakloc_negmeddevs, if any
b=b+(peakloc_negmeddevs-1);
%Calculate cut-off
cutoff = xax(round(b));

%Apply the cutoff
speaks = apeaks(fzapeaks>cutoff);
%Keep mean / sd of resulting distribution
selfzpeaksmean = mean(fzapeaks(fzapeaks>cutoff));
stdfzpeaksmean = std(fzapeaks(fzapeaks>cutoff));

if visualize
    figure
    plot(posmeddevs)
    hold on
    plot(negmeddevs,'k')
    title 'Proportion of positive and negative extremes as a function of threshold'
end



%--------------------------------------------------------------------------
%Step 8: identify IBIs that are too short and remove them

sIBIs = find((diff(speaks)-mean(diff(speaks)))./std(diff(speaks))< -IBIsdco);
remIBI = [];
for i=1:numel(sIBIs)
    csIBI=sIBIs(i);
    
    %Remove either csIBI or csIBI+1
    if sfpulsedat(speaks(csIBI)) > sfpulsedat(speaks(csIBI+1))
        remIBI(i) = csIBI+1;
    else
        remIBI(i) = csIBI;
    end
end

keepIBIbool=ones(size(speaks));
keepIBIbool(remIBI)=0;
speaks=speaks(keepIBIbool==1);


%--------------------------------------------------------------------------
%Step 9: identify IBIs that are too long and try to fix them
lIBIs = find((diff(speaks)-mean(diff(speaks)))./std(diff(speaks))> IBIsdco);

naddPeaks = 0;
addPeaks=[];
for i=1:numel(lIBIs)
    
    stwin = speaks(lIBIs(i));
    enwin = speaks(lIBIs(i)+1);
    
    %How many pulses would normally be in this window?
    nPulses = round((enwin-stwin)/mean(diff(speaks)))-1;

    %Check what peaks there are in that window and what they look like
    canPulses = apeaks(find(apeaks>=stwin&apeaks<=enwin));
    canfz = fzapeaks(find(apeaks>=stwin&apeaks<=enwin));
    
    %If there are any pulses
    if numel(canPulses)>2

        %Where are pulses expected
        for j=1:nPulses
            
            %Location of expected pulse
            cExPulse = round((((canPulses(end)-canPulses(1))/(nPulses+1))*j)+canPulses(1));
            
            %Is there any peak in the vicinity?
            sdcanPulses = find(abs(canPulses-cExPulse)<nSDrep*(std(diff(speaks))));
            
            %Is the highest of those peaks of acceptable shape?
            if numel(sdcanPulses)>0
               
                [maxcanPulse,hsdcanPulses]=max(canfz(sdcanPulses));
                maxcanPulse=maxcanPulse(1); %just to be sure it's only one
                hsdcanPulses=hsdcanPulses(1);
                
                if abs((maxcanPulse-selfzpeaksmean)./stdfzpeaksmean) < nSDtrue
                    
                    %Remember this peak
                    naddPeaks=naddPeaks+1;
                    addPeaks=[addPeaks,canPulses(sdcanPulses(hsdcanPulses))];
                    canPulses(sdcanPulses(hsdcanPulses));
                end            
            end
        end %Loop over expected pulses
    end
end
%Finally add the peaks to the current set
speaks = sortrows(vertcat(speaks,addPeaks'));


%--------------------------------------------------------------------------
%Step 10: Calculate the IBI timeseries
ibits = zeros(numel(sfpulsedat),1);
dspeaks=diff(speaks);

%loop over selected peaks
for i=1:numel(speaks)-1
    ibits(speaks(i)+1:speaks(i+1))= dspeaks(i);
end
%Extrapolate the start and end
ibits(1:find(ibits>0,1)) = ibits(find(ibits>0,1));
ibits=flipud(ibits);
ibits(1:find(ibits>0,1)) = ibits(find(ibits>0,1));
ibits=flipud(ibits);

if visualize
    %Visualize
    figure
    cla
    plot((sfpulsedat.*100)+(mean(ibits)*(1000/sR)),'k')

    hold on
    for i=1:numel(cpeaks)
        plot([cpeaks(i),cpeaks(i)],[600,1400],'g');
    end
    for i=1:numel(speaks)
        plot([speaks(i),speaks(i)],[600,700],'b');
    end
    
    %Plot the IBI timeseries
    plot(ibits.*(1000/sR),'r');
    hold off
    title 'IBI timeseries'
end















