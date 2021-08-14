%% Load and plot experimental data

close all

global EXPDATA

% Load data from MPS experiments. _corrected refers to corrections in 
% measured SD, as described in the manuscript.

fd=cd;
cd DATA

% MPS data
load EXPDATA_G_11mM
load EXPDATA_I_11mM
load EXPDATA_G_11mM_onlyliver

% Human GTT from the study by Dalla Man et al.(2017):
% Man, C. D., Rizza, R. A. and Cobelli, C. (2007). Meal Simulation Model of the Glucose-Insulin
% System. IEEE Trans Biomed Eng 54, 1740–1749. https://doi.org/10.1109/TBME.2007.893506.

load G_mean_Dallaman
load G_SD_Dallaman
load I_mean_Dallaman
load I_SD_Dallaman

% Corrected data with SD

load EXPDATA_G_11mM_corrected
load EXPDATA_I_11mM_corrected
load EXPDATA_G_11mM_onlyliver_corrected

cd(fd)

%Convert input data from table to matrix 
EXPDATA_G=double(table2array(EXPDATA_G_11mM_corrected));
EXPDATA_I=double(table2array(EXPDATA_I_11mM));
EXPDATA_G_onlyliver=double(table2array(EXPDATA_G_11mM_onlyliver_corrected));

G_mean_Dallaman=double(table2array(G_mean_Dallaman));
G_SD_Dallaman=double(table2array(G_SD_Dallaman));
I_mean_Dallaman=double(table2array(I_mean_Dallaman));
I_SD_Dallaman=double(table2array(I_SD_Dallaman));

%Convert time to minutes
G_mean_Dallaman(:,1)=G_mean_Dallaman(:,1)/60;
G_mean_Dallaman(:,1)=G_mean_Dallaman(:,1)-G_mean_Dallaman(1,1);
G_SD_Dallaman(:,1)=G_SD_Dallaman(:,1)/60;
I_mean_Dallaman(:,1)=I_mean_Dallaman(:,1)/60;
I_mean_Dallaman(:,1)=I_mean_Dallaman(:,1)-I_mean_Dallaman(1,1);
I_SD_Dallaman(:,1)=I_SD_Dallaman(:,1)/60;

%Convert glucose concentration to mM
G_mean_Dallaman(:,2)=G_mean_Dallaman(:,2)/18;
G_SD_Dallaman(:,2)=G_SD_Dallaman(:,2)/18;

I_mean_Dallaman(:,2)=I_mean_Dallaman(:,2);
I_SD_Dallaman(:,2)=I_SD_Dallaman(:,2);

%% Data for computational modelling 

% Glucose data from liver+islets
EXPDATA=[];
EXPDATA.time{1} = EXPDATA_G([1 2 4 5 6 7 8 9],1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G([1 2 4 5 6 7 8 9],2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G([1 2 4 5 6 7 8 9],4);      % Time (absolute) (hours)

% Insulin data from liver+islets
EXPDATA.time{2} = EXPDATA_I(1:end,1);    % Time (absolute) (hours)
EXPDATA.mean{2} = EXPDATA_I(1:end,2);    % Time (absolute) (hours)
EXPDATA.SD{2} = EXPDATA_I(1:end,4);      % Time (absolute) (hours)

% Glucose data from only liver
EXPDATA.time{3} = EXPDATA_G_onlyliver([1 3 4 5 7 8 9],1);    % Time (absolute) (hours)
EXPDATA.mean{3} = EXPDATA_G_onlyliver([1 3 4 5 7 8 9],2);    % Time (absolute) (hours)
EXPDATA.SD{3} = EXPDATA_G_onlyliver([1 3 4 5 7 8 9],4);      % Time (absolute) (hours)

% Compute SEM for data in Dalla Man et al.(2017)

SD_G=G_SD_Dallaman(:,2)-G_mean_Dallaman(:,2);
SEM_G=SD_G/sqrt(204);

SD_I=I_SD_Dallaman(:,2)-I_mean_Dallaman(:,2);
SEM_I=SD_I/sqrt(204);




