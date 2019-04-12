%----------------------------------------------------------------------
%BrainAmpConverter v1.2 defaults
%EJH 2011-15
%----------------------------------------------------------------------

%Calibration setting: How many microVolts (the measured unit) corresponds
%to 1 microSiemens? When calibration has not been done online, this value
%usually has to be 25. If it is already done, set this to 1.
uVperuSiem = 1;              

%Canonical pulse window width: Change this setting to determine the width
%(in s) of the canonical pulse wave that is used for classification of
%peaks. Mostly higher values (~8) give more stable results because that way
%the context around a peak is taken into account. Use lower settings (eg 1
%s) for better quality data.
canpulsewin = 1;             %canonical pulse window width (s)

%Enter varieties of response codes for detecting responses (and scan
%triggers). If brainampconverter cannot find responses/triggers which you
%know to be present, look into the brainampfile to see what code identifies
%them. Both 'A' and 'Response' have been used at DCCN.
ResponseCodes = {'A','Response','F'};  %Triggers and button press responses have on of these evt.type code. Add more to work more general.
%Same setting for detecting stimulus codes (when a stimulus is presented).
StimulusCodes = {'Stimulus','S','C'};   %Stimulus presentations have this evt.type code
avoidOverwrite=true;         %Warning when overwriting puls/mat files (false/true)
HRdsFreq = 100;              %Downsample heart rate to 100 Hz
SCRdsFreq = 200;             %Downsample SCR to 200 Hz (because that is default for Autonomate)
extraWindow = 10;            %Extra time window around fMRI runs (for Retroicor)
TriggerBit = 1;              %Trigger is in bit 1 (therefore odd response numbers are triggers)
visualize = false;           %Visualization on/off (for debugging)
SCR_MAFwidth = .25;          %Moving average filter width for SCR
FullFileConversion = 1;      %Convert full file (1) or not (0)

%Channel names to identify the correct recordings (listing all possible names)                             
SCRchannelnames = {'GSR','SCR','EDA','skinconductance'}; %Name of SCR channel in brainamp file (add any option)
HRchannelnames = {'HR','pulse','PPU'};                   %Name of SCR channel in brainamp file (add any option)
Respchannelnames = {'Resp','Respiration','Breath'};               %Name of SCR channel in brainamp file (add any option)
