%Sasha Kramer
%20200320
%UCSB IGPMS

%%%Matrix-based community detection with HPLC pigments:

%*it is highly recommended that you check out the analyses/QC in HPLCcluster 
%and HPLC_EOF before continuing with this code!*%
%https://github.com/sashajane19/HPLCcluster 
%https://github.com/sashajane19/HPLC_EOF 

%Map the directory where you will load your data:
cd /Users/skramer/Documents/UCSB/Research/Data/NAAMES/HPLC

%Load your samples (formatted here as a .mat file):
load N1234_HPLC_surf_nodups.mat

%Remove unwanted pigments (here degradation pigments were already removed 
%[see: HPLCcluster]) so we are removing Tchlb,Tchlc,ABcaro,MVchla,Lutein)
deg =[2,3,4,13,19];
NAAMES_HPLCd = NAAMES_HPLCall;
NAAMES_HPLCd(:,deg) = [];
clear deg

%Normalize all pigment concentrations to total chlorophyll-a
normchl = NAAMES_HPLCd(:,2:end)./NAAMES_HPLCd(:,1); 

%Calculate a correlation matrix for all samples:
corrpig_chlN = corrcoef(normchl'); 

%Apply Weighted Gene Co-Expression Network Analysis (WCGNA; Zhang and
%Horvath, 2005):
abs_CPn = abs(corrpig_chlN); 
WCNA_chlN = abs_CPn.^6; %

for i = 1:length(normchl)
    WCNA_chlN(i,i) = 0;
end
clear i

%GenLouvain method (Jeub et al., 2011-2019):
%Find here: http://netwiki.amath.unc.edu/GenLouvain/GenLouvain 
%Add path to where you saved the files:
addpath(genpath('/Users/skramer/Documents/MATLAB/GenLouvain-master'))
[B,twom] = modularity(WCNA_chlN);
[S1,Q] = genlouvain(B);
mod1 = Q/twom;
clear B twom Q

%Brain connectivity toolbox method (Rubinov and Sporns, 2010):
%Find here: https://sites.google.com/site/bctnet/Home
%Add path to where you saved the files: 
addpath(genpath('/Users/skramer/Documents/MATLAB/connectivity_toolbox'))
[S2,mod2] = modularity_und(WCNA_chlN);

%There are other functions with both of these packages that can help to
%visualize the data, including correlation matrices of samples sorted by
%community results