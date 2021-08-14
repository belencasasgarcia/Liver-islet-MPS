%% Load and plot experimental data

N=10;           %Number of biological replicates in the experiment
L_GTT=48;       %GTT length (hours)

global EXPDATA

fd=cd;
cd DATA

load EXPDATA_G_11mM
load EXPDATA_I_11mM
load EXPDATA_G_5p5mM
load EXPDATA_I_5p5mM
load EXPDATA_G_2p8mM
load EXPDATA_I_2p8mM

cd(fd)

%Convert input data from table to matrix 
EXPDATA_G=double(table2array(EXPDATA_G_11mM));
EXPDATA_I=double(table2array(EXPDATA_I));
EXPDATA_G_5p5mM=double(table2array(EXPDATA_G_5p5mM));
EXPDATA_I_5p5mM=double(table2array(EXPDATA_I_5p5mM));
EXPDATA_G_2p8mM=double(table2array(EXPDATA_G_2p8mM));
EXPDATA_I_2p8mM=double(table2array(EXPDATA_I_2p8mM));

%% Data for modelling 

% Glucose data in the liver-islet co-culture under hyperglycemia (11 mM)
EXPDATA=[];
EXPDATA.time{1} = EXPDATA_G([2 3 4 5 6 8 9],1);    % Time (absolute) (hours)
EXPDATA.mean{1} = EXPDATA_G([2 3 4 5 6 8 9],2);    % Time (absolute) (hours)
EXPDATA.SD{1} = EXPDATA_G([2 3 4 5 6 8 9],4);      % Time (absolute) (hours)

% Insulin data from liver+islets
EXPDATA.time{2} = EXPDATA_I([2 3 4 6 7 8 9],1);    % Time (absolute) (hours)
EXPDATA.mean{2} = EXPDATA_I([2 3 4 6 7 8 9],2);    % Time (absolute) (hours)
EXPDATA.SD{2} = EXPDATA_I([2 3 4 6 7 8 9],4);      % Time (absolute) (hours)

% Glucose data in the liver-islet co-culture under normoglycemia (5.5 mM)

% Glucose data 5.5 mM 
EXPDATA.time{5} = EXPDATA_G_5p5mM([2 4 5],1);    % Time (absolute) (hours)
EXPDATA.mean{5} = EXPDATA_G_5p5mM([2 4 5],2);    % Time (absolute) (hours)
EXPDATA.SD{5} = EXPDATA_G_5p5mM([2 4 5],4);      % Time (absolute) (hours)

% Insulin data 5.5 mM 
EXPDATA.time{6} = EXPDATA_I_5p5mM([2 4 5],1);    % Time (absolute) (hours)
EXPDATA.mean{6} = EXPDATA_I_5p5mM([2 4 5],2);    % Time (absolute) (hours)
EXPDATA.SD{6} = EXPDATA_I_5p5mM([2 4 5],4);      % Time (absolute) (hours)

% Glucose data in the liver-islet co-culture under hypoglycemia (2.8 mM)

% Glucose data 2.8 mM 
EXPDATA.time{7} = EXPDATA_G_2p8mM([1 3 4 5],1);    % Time (absolute) (hours)
EXPDATA.mean{7} = EXPDATA_G_2p8mM([1 3 4 5],2);    % Time (absolute) (hours)
EXPDATA.SD{7} = EXPDATA_G_2p8mM([1 3 4 5],4);      % Time (absolute) (hours)

% Insulin data 2.8 mM 
EXPDATA.time{8} = EXPDATA_I_2p8mM([3 4 5],1);    % Time (absolute) (hours)
EXPDATA.mean{8} = EXPDATA_I_2p8mM([3 4 5],2);    % Time (absolute) (hours)
EXPDATA.SD{8} = EXPDATA_I_2p8mM([3 4 5],4);      % Time (absolute) (hours)

plotExperimentalData(EXPDATA)

