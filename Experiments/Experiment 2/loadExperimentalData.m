% Load experimental data

N=10;           %Number of biological replicates in the experiment
L_GTT=48;       %GTT length (hours)

global EXPDATA

% Load data from MPS experiments. _corrected refers to corrections in 
% measured SD, as described in the manuscript.

fd=cd;
cd DATA

load EXPDATA_G_11mM
load EXPDATA_I_11mM
load EXPDATA_G_11mM_onlyliver
load EXPDATA_I_11mM_onlyislets

load EXPDATA_G_11mM_corrected
%load EXPDATA_I_11mM_corrected
load EXPDATA_G_11mM_onlyliver_corrected

% Glucose and insulin measurements in the individual compartments

load G_CC_liver_overtime
load G_CC_islets_overtime
load I_CC_liver_overtime
load I_CC_islets_overtime

load G_CC_liver_overtime_corrected
load G_CC_islets_overtime_corrected
load I_CC_liver_overtime_corrected

load G_raw

cd(fd)

%Convert input data from table to matrix 

EXPDATA_G=double((EXPDATA_G_11mM_corrected));
EXPDATA_I=double((EXPDATA_I_11mM));
EXPDATA_G_onlyliver=double((EXPDATA_G_11mM_onlyliver_corrected));
EXPDATA_G_liver_overtime=double(table2array(G_CC_liver_overtime));
EXPDATA_G_islets_overtime=double(table2array(G_CC_islets_overtime));
EXPDATA_I_liver_overtime=double(table2array(I_CC_liver_overtime));
EXPDATA_I_islets_overtime=double(table2array(I_CC_islets_overtime));

EXPDATA_G_liver_overtime_corrected=double(table2array(G_CC_liver_overtime_corrected));
EXPDATA_G_islets_overtime_corrected=double(table2array(G_CC_islets_overtime_corrected));
EXPDATA_I_liver_overtime_corrected=double(table2array(I_CC_liver_overtime_corrected));

%% Data for modelling 

% Estimation using first and last GTT

% Glucose data from liver+islets
EXPDATA=[];
EXPDATA.time{1} = EXPDATA_G([1 2 3 4 9 10 11 12],1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G([1 2 3 4 9 10 11 12],2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G([1 2 3 4 9 10 11 12],4);      % Time (absolute) (hours)

% Insulin data from liver+islets
EXPDATA.time{2} = EXPDATA_I([1 2 3 4 9 10 11],1);    % Time (absolute) (hours)
EXPDATA.mean{2} = EXPDATA_I([1 2 3 4 9 10 11],2);    % Time (absolute) (hours)
EXPDATA.SD{2} = EXPDATA_I([1 2 3 4 9 10 11],4);      % Time (absolute) (hours)

% Glucose data from only liver
EXPDATA.time{3} = EXPDATA_G_onlyliver([1 2 4 9 10 11 12],1);    % Time (absolute) (hours)
EXPDATA.mean{3} = EXPDATA_G_onlyliver([1 2 4 9 10 11 12],2);    % Time (absolute) (hours)
EXPDATA.SD{3} = EXPDATA_G_onlyliver([1 2 4 9 10 11 12],4);      % Time (absolute) (hours)

% Measurements of glucose and insulin over time

% Glucose data in the co-culture measured in the liver compartment over
% time
EXPDATA.time{4} = EXPDATA_G_liver_overtime_corrected(2:end-2,1);
EXPDATA.mean{4} = EXPDATA_G_liver_overtime_corrected(2:end-2,2);
EXPDATA.SD{4} = EXPDATA_G_liver_overtime_corrected(2:end-2,4);

% Glucose data in the co-culture measured in the islets compartment over
% time
EXPDATA.time{5} = EXPDATA_G_islets_overtime_corrected(2:end-2,1);
EXPDATA.mean{5} = EXPDATA_G_islets_overtime_corrected(2:end-2,2);
EXPDATA.SD{5} = EXPDATA_G_islets_overtime_corrected(2:end-2,4);

% Insulin data in the co-culture measured in the liver compartment over
% time
EXPDATA.time{6} = EXPDATA_I_liver_overtime_corrected(2:end-2,1);
EXPDATA.mean{6} = EXPDATA_I_liver_overtime_corrected(2:end-2,2);
EXPDATA.SD{6} = EXPDATA_I_liver_overtime_corrected(2:end-2,4);

% Insulin data in the co-culture measured in the islets compartment over
% time
EXPDATA.time{7} = EXPDATA_I_islets_overtime(2:end-2,1);
EXPDATA.mean{7} = EXPDATA_I_islets_overtime(2:end-2,2);
EXPDATA.SD{7} = EXPDATA_I_islets_overtime(2:end-2,4);

% Raw data for glucose
EXPDATA.G_raw=G_raw;

