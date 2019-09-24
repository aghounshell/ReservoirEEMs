% code written by Rose M. Cory
% code adapted by D. Scott, 20101030
% Updated by A.G. Hounshell, 20190924

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
folder = 'C:\Users\ahounshell\Documents\VT_PostDoc\EEMs\20190703\';
sample = '1';
fld = '\';                  % / for mac; \ for PC
% Calculate Raman Area using the Raman Scan collected on the same day as
% analysis
raman_area = 1.2736e7   %enter respective raman area for each EEM, intensities will be normalized to this value
dilution_factor = 1; %use this if you want
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
temp = strcat(folder,folder_f,fld,sample,sampleext);   % file name
importfile(temp);
A = data;
% A = importdata(temp);
clear colheaders textdata data temp;
ifile = strcat(sample,samplecor);    %whatever you want to name the output file
%
%% Step 3. Read in blank EEM
%
temp = strcat(folder,folder_b,fld,fnm_blank_eem);    
importfile(temp);
Ab = data;
clear colheaders textdata data temp;
%
%% Step 4. Read in Absorbance file for corrections
%
temp = strcat(folder,folder_a,fld,sample,sampleabs);    
%abs = csvread('/Volumes/durelles/science/flourescence/2010/test/absorbance/pony_dil.csv',0,0);   %read in absorbance file for inner filter correction (see below for the actual correction)
abs = csvread(temp,0,0);
temp = strcat(folder,folder_a,fld,blankabs);    
%absb = csvread('/Volumes/durelles/science/flourescence/2010/test/absorbance/blank.csv',0,0);   %read in absorbance blank file for inner filter correction (see below for the actual correction)
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
Asize = size(A);
emissionLen = Asize(1) ;
emission =   A(1:emissionLen, 1);  
% Cut emission wavelengths from A and Ab (blank matrix)
A=A(:,2:44); %cuts em wavelengths from 'A'
Ab = Ab(:,2:44)
% Substract blank EEM from sample EEM
A = A-Ab
clear Ab
Asize = size(A);
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
z = A(1:ylen, 1:xlen); %redefines 'A' as 'z'
zi = interp2(x, y, z, xi, yi, 'spline'); %interpolation
%
%READ IN THE EX AND EM CORRECTION FILES (NOTE: EX FILE IS INVERTED....IF
%YOU DON'T INVERT THEN CHANGE EQUATION BELOW (LINE 78) TO DIVIDE BY EX CORRECTION)
% 
% DOUBLE CHECK THIS!!!!! WAS WRONG ORIGINALLY - NEED TO USE UPDATED
% CORRECT FILES
%
temp = strcat(folder,folder_e,fld,'mcorrect_f4_350_550_1.xls');    
MC = xlsread(temp);
temp = strcat(folder,folder_e,fld,'xcorrect_f4_240_450_12_5.xls');   
XC = xlsread(temp);
%
%APPLYING EX AND EM CORRECTIONS WITH MATRIX ALGEBRA
%
X=diag(XC); %creates ex correction matrix
%
Y=diag(MC); %creates em correcting matrix
%
zi=[[zi*X]'*Y]'; %applies corrections
%
%CALCULATES FI (EX=370; EM470/EM520)
FI=zi(121,53)/zi(171,53); 
%
%NORMALIZES CORRECTED DATA TO RAMAN AREA
%
zir=zi/raman_area; 
%
% ---------------------------------
%INNERFILTER CORRECTION
abs_ex1=abs_ex1(:,2:2); %cuts wavelength labels and stnd. dev. from file
%
ao = 190:5:850; %defines wavelength range
%
ai = 190:2.5:850; %defines wavelength range to interpolate to
%
iabs_ex1 = interp1(ao,abs_ex1,ai); %interpolates to 2.5 nm
%
iabs_ex1 = (iabs_ex1)'; %moves data from columns to rows
%
ex_abs=iabs_ex1(21:105,:); %selects data from 240-450 (ex range)
%
ex_abs=ex_abs'; %moves data from rows to columns
%
abs = flipud(abs);
em_abs=abs(161:361,:); %selects data from 350-550 (em range)
%
em_abs=em_abs(:,2:2); %cuts wavelength labels and stnd. dev. from file
%
for i=1:length(em_abs)
    for j=1:length(ex_abs)
        IFC(i,j)=ex_abs(j)+em_abs(i); %defines 'IFC' as the sum of the ex and em wavelengths for all ex/em pairs 
    end
end
czir=zir.*10.^(0.5*IFC); %applies inner filter correction
%
% CALCULATES FI AFTER ALL CORRECTIONS HAVE BEEN APPLIED
FInew=czir(121,53)/czir(171,53); %calculates FI ex=370 nm em470/em520
%
%SAVE THE RAMAN NORMALIZED AND EX, EM, AND INNER FILTERED CORRECTED EEM MATRIX
%
pathname = strcat(folder,folder_d,fld);
%
for i=1:length(ifile)
    pathname(length(pathname) + 1) = ifile(i);
end%

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

em = 350:1:550; 

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

pathname = strcat(folder,folder_c,fld);

for i=1:length(ifile)

    pathname(length(pathname) + 1) = ifile(i);

end

pathnamelength = length(pathname);

pathname(pathnamelength + 1: pathnamelength + 4) = '.tif';

saveas(gcf, pathname, 'tif');