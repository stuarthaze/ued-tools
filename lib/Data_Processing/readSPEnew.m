% readSPE.m
function [image, header] = readSPEnew(dirPath, filename)
% readSPE.m  Read princeton instrument *.SPE 2.x image files
%
% image, header = readSPE(filename)
%               = readSPE(dirPath, filename)
%
% Inputs:  filename - path and filename (string)
%           dirPath - optional directory path (string)
% 
% Outputs:    image - 3D array of image(s)
%            header - specific header info (struct)
% 
% Image is returned as a 3D array, where each image is the first 2
% dimensions, and successive images are stored along the 3rd dimension.
% Modify the contents of the output header struct to suit your uses.
%
% The image is stored as it is shown in WinVIEW; the first two dimensions 
% are stored as [pixel,stripe]
% 
% Only minimal error checking, and only the uint16 pixel datatype has been tested, so
% be more cautious than usual if using other datatypes.
%
% The SPE 2.x header byte structure is taken from Appendix A of PI's documentation,
%    ftp://ftp.princetoninstruments.com/public/Manuals/Princeton%20Instruments/SPE%203.0%20File%20Format%20Specification.pdf
%  
% This method is based off the following sources:
%  
%  Matt Clarke's Python Script Description:
%   http://www.scipy.org/Cookbook/Reading_SPE_files
%
%  Stuart B. Wilkins's Python Code:
%   https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/ccd/files.py 
%
%  ImageJ's plugin to import SPE files, written in Java:
%   http://rsbweb.nih.gov/ij/
%  

% Author: Carl Hall
% Date: June 2012

% Modification History
%  April 2012 - Original Code
%   June 2012 - Modified to maintain output array as same datatype as SPE
%               file
%  March 2017 - Incorporate comments from Tremayne and Shijie to use fullfile()
%               and fclose()
%             - Update description for source of header byte offsets
%             - Use decimal offsets instead of hex to match documentation
%             - Add ability to output some of the header info

%% SPE 2.x Header byte Structure (not fully checked) (little endian)
% 
% Type    Variable                     Offset   Description 
%  int16   ControllerVersion               0     Hardware Version
%  int16   LogicOutput                     2     Definition of Output BNC
%  uint16  AmpHiCapLowNoise                4     Amp Switching Mode
%  uint16  xDimDet                         6     Detector x dimension of chip.
%  int16   mode                            8     timing mode
%  float   exp_sec                        10     alternative exposure, in sec.
%  int16   VChipXdim                      14     Virtual Chip X dim
%  int16   VChipYdim                      16     Virtual Chip Y dim
%  uint16  yDimDet                        18     y dimension of CCD or detector.
%  int8    date[DATEMAX]                  20     date
%  int16   VirtualChipFlag                30     On/Off
%  int8    Spare_1[2]                     32 
%  int16   noscan                         34     Old number of scans - should always be -1
%  float   DetTemperature                 36     Detector Temperature Set
%  int16   DetType                        40     CCD/DiodeArray type
%  uint16  xdim                           42     actual # of pixels on x axis
%  int16   stdiode                        44     trigger diode
%  float   DelayTime                      46     Used with Async Mode
%  uint16  ShutterControl                 50     Normal, Disabled Open, Disabled Closed
%  int16   AbsorbLive                     52     On/Off
%  uint16  AbsorbMode                     54     Reference Strip or File
%  int16   CanDoVirtualChipFlag           56     T/F Cont/Chip able to do Virtual Chip
%  int16   ThresholdMinLive               58     On/Off
%  float   ThresholdMinVal                60     Threshold Minimum Value
%  int16   ThresholdMaxLive               64     On/Off
%  float   ThresholdMaxVal                66     Threshold Maximum Value
%  int16   SpecAutoSpectroMode            70     T/F Spectrograph Used 
%  float   SpecCenterWlNm                 72     Center Wavelength in Nm
%  int16   SpecGlueFlag                   76     T/F File is Glued
%  float   SpecGlueStartWlNm              78     Starting Wavelength in Nm
%  float   SpecGlueEndWlNm                82     Starting Wavelength in Nm
%  float   SpecGlueMinOvrlpNm             86     Minimum Overlap in Nm
%  float   SpecGlueFinalResNm             90     Final Resolution in Nm
%  int16   PulserType                     94     0=None, PG200=1, PTG=2, DG535=3
%  int16   CustomChipFlag                 96     T/F Custom Chip Used
%  int16   XPrePixels                     98     Pre Pixels in X direction
%  int16   XPostPixels                   100    Post Pixels in X direction
%  int16   YPrePixels                    102     Pre Pixels in Y direction
%  int16   YPostPixels                   104     Post Pixels in Y direction
%  int16   asynen                        106     asynchron enable flag  0 = off
%  int16   datatype                      108     pixel datatype: 0 = float (4 bytes)
%                                                                1 = int32 (4 bytes)
%                                                                2 = int16 (2 bytes)
%                                                                3 = uint16 (2 bytes)
%                                                                8 = uint32 (4 bytes) 
%  int16   PulserMode                    110     Repetitive/Sequential
%  uint16  PulserOnChipAccums            112     Num PTG On-Chip Accums
%  uint32  PulserRepeatExp               114     Num Exp Repeats (Pulser SW Accum)
%  float   PulseRepWidth                 118     Width Value for Repetitive pulse (usec)
%  float   PulseRepDelay                 122     Width Value for Repetitive pulse (usec)
%  float   PulseSeqStartWidth            126     Start Width for Sequential pulse (usec)
%  float   PulseSeqEndWidth              130     End Width for Sequential pulse (usec)
%  float   PulseSeqStartDelay            134     Start Delay for Sequential pulse (usec)
%  float   PulseSeqEndDelay              138     End Delay for Sequential pulse (usec)
%  int16   PulseSeqIncMode               142     Increments: 1=Fixed, 2=Exponential
%  int16   PImaxUsed                     144     PI-Max type controller flag
%  int16   PImaxMode                     146     PI-Max mode
%  int16   PImaxGain                     148     PI-Max Gain
%  int16   BackGrndApplied               150     1 if background subtraction done
%  int16   PImax2nsBrdUsed               152     T/F PI-Max 2ns Board Used
%  uint16  minblk                        154     min. # of strips per skips
%  uint16  numminblk                     156     # of min-blocks before geo skps
%  int16   SpecMirrorLocation[2]         158     Spectro Mirror Location, 0=Not Present
%  int16   SpecSlitLocation[4]           162     Spectro Slit Location, 0=Not Present 
%  int16   CustomTimingFlag              170     T/F Custom Timing Used
%  int8    ExperimentTimeLocal[TIME MAX] 172     Experiment Local Time as hhmmss\0 
%  int8    ExperimentTimeUTC[TIME MAX]   179     Experiment UTC Time as hhmmss\0 
%  int16   ExposUnits                    186     User Units for Exposure
%  uint16  ADCoffset                     188     ADC offset
%  uint16  ADCrate                       190     ADC rate
%  uint16  ADCtype                       192     ADC type
%  uint16  ADCresolution                 194     ADC resolution
%  uint16  ADCbitAdjust                  196     ADC bit adjust
%  uint16  gain                          198     gain
%  int8    Comments[5][COMMENTMAX]       200     File Comments 
%  uint16  geometric                     600     geometric ops: rotate 0x01,reverse 0x02, flip 0x04 
%  int8    xlabel[LABELMAX]              602     intensity display string
%  uint16  cleans                        618     cleans
%  uint16  NumSkpPerCln                  620     number of skips per clean.
%  int16   SpecMirrorPos[2]              622     Spectrograph Mirror Positions
%  float   SpecSlitPos[4]                626     Spectrograph Slit Positions
%  int16   AutoCleansActive              642     T/F
%  int16   UseContCleansInst             644     T/F
%  int16   AbsorbStripNum                646     Absorbance Strip Number
%  int16   SpecSlitPosUnits              648     Spectrograph Slit Position Units
%  float   SpecGrooves                   650     Spectrograph Grating Grooves
%  int16   srccmp                        654     number of source comp.diodes
%  uint16  ydim                          656     y dimension of raw data.
%  int16   scramble                      658     0=scrambled,1=unscrambled
%  int16   ContinuousCleansFlag          660     T/F Continuous Cleans Timing Option
%  int16   ExternalTriggerFlag           662     T/F External Trigger Timing Option
%  int32   lnoscan                       664     Number of scans (Early WinX)
%  int32   lavgexp                       668     Number of Accumulations
%  float   ReadoutTime                   672     Experiment readout time
%  int16   TriggeredModeFlag             676     T/F Triggered Timing Option
%  uint64  XML Offset                    678     Starting location of the XML footer
%  int8    sw_version[FILEVERMAX]        688     Version of SW creating this file 
%  int16   type                          704     1 = new120 (Type II)  
%                                                2 = old120 (Type I) 
%                                                3 = ST130
%                                                4 = ST121
%                                                5 = ST138
%                                                6 = DC131 (PentaMax)
%                                                7 = ST133 (MicroMax/SpectroMax)
%                                                8 = ST135 (GPIB) 
%                                                9 = VICCD
%                                               10 = ST116 (GPIB)
%                                               11 = OMA3 (GPIB)
%                                               12 = OMA4 
%  int16   flatFieldApplied              706     1 if flat field was applied
%  int8    Spare_3[16]                   708 
%  int16   kin_trig_mode                 724     Kinetics Trigger Mode
%  int8    dlabel[LABELMAX]              726     Data label.
%  int8    Spare_4[436]                  742 
%  int8    PulseFileName[HDRNAMEMAX]    1178     Name of Pulser File with Pulse Widths/Delays (for Z-Slice)
%  int8    AbsorbFileName[HDRNAMEMAX]   1298     Name of Absorbance File (if File Mode)
%  uint32  NumExpRepeats                1418     Number of Times experiment repeated
%  uint32  NumExpAccums                 1422     Number of Time experiment accumulated
%  int16   YT_Flag                      1426     Set to 1 if this file contains YT data
%  float   clkspd_us                    1428     Vert Clock Speed in micro-sec
%  int16   HWaccumFlag                  1432     set to 1 if accum done by Hardware.
%  int16   StoreSync                    1434     set to 1 if store sync used
%  int16   BlemishApplied               1436     set to 1 if blemish removal applied
%  int16   CosmicApplied                1438     set to 1 if cosmic ray removal applied
%  int16   CosmicType                   1440     if cosmic ray applied, this is type
%  float   CosmicThreshold              1442     Threshold of cosmic ray removal.
%  int32   NumFrames                    1446     number of frames in file.
%  float   MaxIntensity                 1450     max intensity of data (future)
%  float   MinIntensity                 1454     min intensity of data future)
%  int8    ylabel[LABELMAX]             1458     y axis label.
%  uint16  ShutterType                  1474     shutter type.
%  float   shutterComp                  1476     shutter compensation time.
%  uint16  readoutMode                  1480     readout mode, full,kinetics, etc
%  uint16  WindowSize                   1482     window size for kinetics only.
%  uint16  clkspd                       1484     clock speed for kinetics & frame transfer
%  uint16  interface_type               1486     computer interface (isa-taxi, pci, eisa, etc.)
%  int16   NumROIsInExperiment          1488     May be more than the 10 allowed in this header (if 0, assume 1)
%  int8    Spare_5[16]                  1490 
%  uint16  controllerNum                1506     if multiple controller system will have controller number data camefrom. This is a future item. 
%  uint16  SWmade                       1508     Which software package created this file 
%  int16   NumROI                       1510     number of ROIs used. if 0 assume 1.
%          ROIinfo[6], contains:
%  uint16   startx                               left x start value.
%  uint16   endx                                 right x value.
%  uint16   groupx                               amount x is binned/grouped in hw.
%  uint16   starty                               top y start value.
%  uint16   endy                                 bottom y value.
%  uint16   groupy                               amount y is binned/grouped in hw.
%                                       1512     ROI 1
%                                       1524     ROI 2
%                                       1536     ROI 3
%                                       1548     ROI 4 
%                                       1560     ROI 5
%                                       1572     ROI 6
%                                       1584     ROI 7
%                                       1596     ROI 8
%                                       1608     ROI 9
%                                       1620     ROI 10
%  char    FlatField[120]               1632     Flat field file name.
%  char    background[120]              1752     background sub. file name.
%  char    blemish[120]                 1872     blemish file name.
%  float   file_header_ver              1992     version of this file header (3.0)
%  int8    YT_Info[1000]                1996     Reserved for YT information 
%  int32   WinView_id                   2996     == 0x01234567L if file created by WinX 
%          Image data                   4100

%% Start of Code

% parse optional input
if nargin>1
    filename = fullfile(dirPath,filename);
else
    filename = dirPath;
end

% Open the file
fd = fopen(filename,'r');
if(fd < 0)
    error('Could not open file, bad filename');
end

% Get the image dimensions and pixel datatype
header = struct('stripDim',getData(fd,   42, 'uint16'),...      %first dim
                'pixelDim',getData(fd,  656, 'uint16'),...      %second dim
                'nDim'    ,getData(fd, 1446, 'uint32'),...      %third dim
                'dataType',getData(fd,  108, 'uint16'));

% Get additional useful data from header
header.exposure = getData(fd, 118, 'float') * 1e3;  %ICCD exposure, in ns
header.MCPGain  = getData(fd, 148, 'int16');        %ICCD gain

% Get the image
fseek(fd, 4100, 'bof');
image = zeros([header.pixelDim,header.stripDim,header.nDim]);
switch header.dataType
    case 0     % single precision float (4 bytes)
        image = single(image);      %maintain datatype in function output
        for i=1:header.nDim
            image(:,:,i) = fread(fd, [header.stripDim,header.pixelDim], 'float32')';
        end
    
    case 1    % long int (4 bytes)
        image = int32(image);
        for i=1:header.nDim
            image(:,:,i) = fread(fd, [header.stripDim,header.pixelDim], 'int32')';
        end
    
    case 2    % short int (2 bytes)
        image = int16(image);
        for i=1:header.nDim
            image(:,:,i) = fread(fd, [header.stripDim,header.pixelDim], 'int16')';
        end
    
    case 3    % short unsigned int (2 bytes)
        image = uint16(image);
        for i=1:header.nDim
            image(:,:,i) = fread(fd, [header.stripDim,header.pixelDim], 'uint16')';
        end

    otherwise
        error('Unknown pixel datatype');
end

fclose(fd);
end

%% 
% getData() reads one piece of data at a specific location
% 
function data = getData(fd, decLoc, dataType)
% Inputs: fd - int    - file descriptor
%     decLoc - int    - location of data relative to beginning of file, decimal
%   dataType - string - type of data to be read
%
fseek(fd, decLoc, 'bof');
data = fread(fd, 1, dataType);
end