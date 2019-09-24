% Script to manipulate Ex and Em Correction files from FluoroMax4
% Correction files obtained as .spc files from instrument
% A Hounshell, 24Sep2019

clear

% Folder: C:\Users\ahoun\Dropbox\Reservoir_EEMs\Instrument_CorrFiles
% Load in .spc files
mcorr = tgspcread('MCorrect.spc');
xcorr = tgspcread('XCorrect.spc');

% Export as .csv files in original range/wavelengths
mcorr.data(:,1) = mcorr.X;
mcorr.data(:,2) = mcorr.Y;
csvwrite('mcorr.csv',mcorr.data);

xcorr.data(:,1) = xcorr.X;
xcorr.data(:,2) = xcorr.Y;
csvwrite('xcorr.csv',xcorr.data);

% Manipulate data to constrain to Single_EEM_correction file
% Xcorr: 240-450 every 2.5 nm
% Interpolate every 0.5 nm
xcorr.half = 240:0.5:450;
xcorr.half = xcorr.half';
xcorr.half(:,2) = interp1(xcorr.data(:,1),xcorr.data(:,2),xcorr.half(:,1));
% Select every 2.5 nm wavelengths from 240-450 nm
icount = 0;
for i = 1:5:421
    icount = icount + 1;
    xcorr.nm2(icount,1) = xcorr.half(i,1);
    xcorr.nm2(icount,2) = xcorr.half(i,2);
end
% Export file
xlswrite('xcorrect_f4_240_450_12_5_corr.xls',xcorr.nm2);

% MCorr: 300-600 every 1 nm
% Interpolate every 1 nm
mcorr.one = 300:1:600;
mcorr.one = mcorr.one';
mcorr.one(:,2) = interp1(mcorr.data(:,1),mcorr.data(:,2),mcorr.one(:,1));
% Export file
xlswrite('mcorrect_f4_300_600_1_corr.xls',mcorr.one);

% Manipulate data to constrain to wavelengths analyzed for AGH EEMs
% XCorr: 240-450 every 5 nm
icount = 0;
for i = 1:5:211
    icount = icount + 1;
    xcorr.nm5(icount,1) = xcorr.data(i,1);
    xcorr.nm5(icount,2) = xcorr.data(i,2);
end
% Export file
xlswrite('xcorr_240_450_5.xls',xcorr.nm5);

% MCorr: 300-600 every 2 nm
% Select data for every 2 nm from interpolated data
icount = 0;
for i = 1:2:301
    icount = icount + 1;
    mcorr.nm2(icount,1) = mcorr.one(i,1);
    mcorr.nm2(icount,2) = mcorr.one(i,2);
end
% Export file
xlswrite('xcorr_300_600_2.xls',mcorr.nm2);

% Data files saved in
% 'C:\Users\ahoun\Dropbox\Reservoir_EEMs\Instrument_CorrFiles'