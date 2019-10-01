% code written by Rose M. Cory
% code adapted by D. Scott, 20101030
% Updated by A.G. Hounshell, 20190924
%   Calcuate Raman_Area from daily Raman Scan
%   Updated Instrument correction files
%   Use 'for loop' to calculate instrument Ex and Em corrections

% Purpose of this code is to correct and plot an EEM. corrected EEMs and
% Figures are saved to set folders. 
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% STEP 1. FOLDER and SAMPLE Names & Raman Area
%
% Note: forward slash on Apple/backward slash on PC
%
% Select folder where data is located
folder = 'C:\Users\ahoun\Dropbox\Reservoir_EEMs\30Sep19';
% Input sample number
sample = '2';
fld = '\';                  % / for mac; \ for PC
% Calculate Raman Area using the Raman Scan collected on the same day as
% analysis
% Select Raman File
[datafile_name, directory_name] = uigetfile({'*.csv'},'Choose file for processing');
cd(directory_name);
data = importdata(datafile_name);
raman.raw = data.data; % save the original file.
% Apply instrument corrections to Raman file
raman.corrfile = xlsread('mcorrect_raman.xls');
for i = 1:86
    raman.corr(i,1) = raman.corrfile(i,2)*raman.raw(i,4)/4.61251545;
end
% Calculate area under Raman peak
raman_area = trapz(smooth(raman.corr(6:64,1)));
% Input dilution factor (use if sample was diluted!)
dilution_factor = 1;
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Folder names & standard files below; do not modify
folder_a = 'absorbance';
folder_b = 'blank_EEM';
folder_c = 'corrected_EEM_figs';
folder_d = 'corrected_scans';;
folder_e = 'instrum_correction'
folder_f = 'raw_EEM';
sampleext = '_eem.csv';
samplecor = '_eemc.csv';
fnm_blank_eem = 'b_eem.csv'
sampleabs = '_abs.csv';
blankabs = 'b_abs.csv';
%
%% STEP 2. READ IN EEM
%
temp = strcat(folder,fld,folder_f,fld,sample,sampleext);   % file name
A = importdata(temp);
clear colheaders textdata data temp;
ifile = strcat(sample,samplecor);    %whatever you want to name the output file
%
%% Step 3. Read in blank EEM
%
temp = strcat(folder,fld,folder_b,fld,fnm_blank_eem); 
Ab = importdata(temp);
clear colheaders textdata data temp;
%
%% Step 4. Read in Absorbance file for corrections
%
temp = strcat(folder,fld,folder_a,fld,sample,sampleabs);    
abs = csvread(temp,0,0);
temp = strcat(folder,fld,folder_a,fld,blankabs);    
absb = csvread(temp,0,0);
abs(:,2) = abs(:,2) - absb(:,2); % subtract blank absorbance
clear absb
abs_ex1 = abs(1:5:end, :); %takes abs values every 5nm
abs_ex1 = flipud(abs_ex1); % flips 
%
%% Step 5. Define excitation and emission wavelengths. Ensure correct!!
%
% Define excitation wavelengths
%
excitation = (240:5:450); % defines ex wavelengths
%
% Define emission wavelengths, cut labels from matrix A (uncorrected EEM)
% Cut emission wavelengths from A and Ab (blank matrix)
A.data = A.data(2:152,:); %cuts em wavelengths from 'A'
Ab.data = Ab.data(2:152,:);
Asize = size(A.data);
emissionLen = Asize(1);
emission =   300:2:600;  % May need to change based on wavelengths used
% Substract blank EEM from sample EEM
A.data = A.data-Ab.data;
clear Ab
%
%REDEFINES EX AND EM AS X AND Y
%
ylen = Asize(1);
xlen = Asize(2);
y = emission;
x = excitation;
xend = x(xlen);
yend = y(ylen);
%
%INTERPOLATING THE DATA
%
[xi, yi] = meshgrid(x(1):2.5:xend, y(1):1:yend); %defines wavelengths to interpolate to (ex by 2.5 and em by 1)
z = A.data(1:ylen, 1:xlen); %redefines 'A' as 'z'
zi = interp2(x, y, z, xi, yi, 'spline'); %interpolation
%
%READ IN THE EX AND EM CORRECTION FILES
%
temp = strcat(folder,fld,folder_e,fld,'mcorrect_f4_300_600_1_corr.xls');    
MC = xlsread(temp);
temp = strcat(folder,fld,folder_e,fld,'xcorrect_f4_240_450_12_5_corr.xls');   
XC = xlsread(temp);
%
%APPLYING EX AND EM CORRECTIONS: Using Matlab for loop
%
for i = 1:301
    for j = 1:85
       zi(i,j) = zi(i,j)/XC(j,2)*MC(i,2); 
    end
end
%
%CALCULATES FI (EX=370; EM470/EM520)
FI=zi(171,53)/zi(221,53); 
%
%NORMALIZES CORRECTED DATA TO RAMAN AREA
%
zir=zi/raman_area; 
%
% ---------------------------------
%INNERFILTER CORRECTION
ao = 190:5:850; %defines wavelength range
%
ai = 190:2.5:850; %defines wavelength range to interpolate to
%
iabs_ex1(:,1) = flipud(ai');
iabs_ex1(:,2) = interp1(ao,abs_ex1(:,2),ai); %interpolates to 2.5 nm
%
ex_abs=iabs_ex1(161:245,:); %selects data from 240-450 (ex range)
ex_abs=flipud(ex_abs);
%
abs = flipud(abs);
em_abs=abs(251:551,:); %selects data from 300-600 (em range)
em_abs=flipud(em_abs);
%
for i=1:length(em_abs)
    for j=1:length(ex_abs)
        IFC(i,j)=ex_abs(j,2)+em_abs(i,2); %defines 'IFC' as the sum of the ex and em wavelengths for all ex/em pairs 
    end
end
czir=zir.*10.^(0.5*IFC); %applies inner filter correction
%
% CALCULATES FI AFTER ALL CORRECTIONS HAVE BEEN APPLIED
FInew=czir(171,53)/czir(221,53); %calculates FI ex=370 nm em470/em520
%
%SAVE THE RAMAN NORMALIZED AND EX, EM, AND INNER FILTERED CORRECTED EEM MATRIX
%
pathname = strcat(folder,fld,folder_d,fld);
%
for i=1:length(ifile)
    pathname(length(pathname) + 1) = ifile(i);
end
%
pathnamelength = length(pathname);
%
pathname(pathnamelength + 1: pathnamelength + 4) = '.xls';
%
save(pathname, 'czir', '-ascii', '-double', '-tabs');
%
%
%
%% PLOTS THE EEM 3D CONTOUR STYLE, 
ex = 240:2.5:450;

em = 300:1:600; 

%A=czir';
A=czir;
A=A*dilution_factor
A=A'
%
ex = 240:2.5:400;
A = A(1:65,:)
figure; %to see the figure, run the restof the code all at once from this point on 

%contourf(em,ex,A,500);  
contourf(em,ex,A,300); % change value of 500 to obtain resolution in plot
handle = gca;

set(handle,'fontsize', 14);

%colormap(cool);
%colormap(jet);
colormap(copper);

caxis([.001, 1.2]); %change the scale of the intensity axis

caxis('manual');
xlabel('Emission (nm)')
ylabel('Excitation (nm)')

H = colorbar('vert');

set(H,'fontsize',14);

pathname = strcat(folder,fld,folder_c,fld);

for i=1:length(ifile)

    pathname(length(pathname) + 1) = ifile(i);

end

pathnamelength = length(pathname);

pathname(pathnamelength + 1: pathnamelength + 4) = '.tif';

saveas(gcf, pathname, 'tif');