% Optimize model parameters using the data EXPDATA

clear
clc
close all

format long
format compact

% Declare global variables

global EXPDATA
global modelName
global icOrig0
global pNames
global param
global OPTIONS
global COSTOPTIONS
global SIMTIME      % Time axis for simulation
global time_highres % Time axis for plotting simulations 
global OPTIME       % Time axis for optimization
global N_VAR        % Number of variables to optimize
global OPTVAR_INDEX % Indexes of variables to be optimized
global CUTOFF
global dataPoints
global FID
global Optparam
global parIndex


%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

modelName ='MPSmodel';

optModel = IQMmodel(strcat(modelName,'.txt')); 
IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

%% Parameter indexes

loadParameterIndexes

%% Optimization settings

% Simulation time
time_end=360;

SIMTIME=[0:0.1:time_end];

time_highres=[0:0.001:time_end];

% Change experimental data based on the simulation time

% Find experimental values corresponding to the simulation time

for i=1:size(EXPDATA.time,2)
    time_indexes=find(EXPDATA.time{i}<=time_end);
    EXPDATA.time{i}=EXPDATA.time{i}(time_indexes);
    EXPDATA.mean{i}=EXPDATA.mean{i}(time_indexes);
    EXPDATA.SD{i}=EXPDATA.SD{i}(time_indexes);
end

plotExperimentalData(EXPDATA)

Optparam=param;  %Optimal parameters

OPTVAR_INDEX=[1 2 3];    %Experimental data to optimize against. 
OPTIONS.lowbounds=[];
OPTIONS.highbounds=[];

OPTIONS.lowbounds(1)=1e-2;                % EG0
OPTIONS.lowbounds(2)=1e-9;                % S_i
OPTIONS.lowbounds(3)=1e6;                 % Sigma
OPTIONS.lowbounds(4)=1e-9;                % CL_i_hep
OPTIONS.lowbounds(5)=1e-9;                % Imax_SI
OPTIONS.lowbounds(6)=100;                 % EC50_SI
OPTIONS.lowbounds(7)=1;                   % kv
OPTIONS.lowbounds(8)=48;                  % Alpha
OPTIONS.lowbounds(9)=-2;                  % delta_G_1_11mM
OPTIONS.lowbounds(10)=-2;                 % delta_G_13_11mM
OPTIONS.lowbounds(11)=0;                  % delta_I_1_11mM
OPTIONS.lowbounds(12)=-2;                 % delta_G_1_liver

OPTIONS.lowbounds=OPTIONS.lowbounds';

% High parameter bounds

OPTIONS.highbounds(1)=1e9;                % EG0
OPTIONS.highbounds(2)=0.1;                % S_i
OPTIONS.highbounds(3)=1e9;                % Sigma
OPTIONS.highbounds(4)=200;                % CL_i_hep
OPTIONS.highbounds(5)=1;                  % Imax_SI
OPTIONS.highbounds(6)=1e5;                % EC50_SI
OPTIONS.highbounds(7)=5;                  % kv
OPTIONS.highbounds(8)=500;                % Alpha
OPTIONS.highbounds(9)=2;                  % delta_G_1_11mM
OPTIONS.highbounds(10)=2;                 % delta_G_13_11mM
OPTIONS.highbounds(11)=500;               % delta_I_1_11mM
OPTIONS.highbounds(12)=2;                 % delta_G_1_liver

OPTIONS.highbounds=OPTIONS.highbounds';

startGuess(1)=1;                          % EG0
startGuess(2)=0.01;                       % S_i
startGuess(3)=4.8e6;                      % Sigma
startGuess(4)=20;                         % CL_i_hep 
startGuess(5)=0.92;                       % Imax_SI  
startGuess(6)=550;                        % EC50_SI 
startGuess(7)=5;                          % kv
startGuess(8)=392;                        % Alpha
startGuess(9)=0;                          % delta_G_1_11mM
startGuess(10)=0;                         % delta_G_13_11mM
startGuess(11)=0;                         % delta_I_1_11mM
startGuess(12)=0;                         % delta_G_1_liver


X=startGuess;
 
OPTIONS.index_Optpar=boolean(sum([parIndex.i_EG0_hep...
    parIndex.i_S_i parIndex.i_Sigma parIndex.i_CL_i_hep...
    parIndex.i_Imax_SI parIndex.i_EC50_SI...
    parIndex.i_kv parIndex.i_Alpha...
    parIndex.i_delta_G_1_11mM parIndex.i_delta_G_13_11mM parIndex.i_delta_I_1_11mM...
    parIndex.i_delta_G_1_liver],2)); %Indexes of the parameters to optimize


%Set parameter values manually based on a-priori knowledge
Optparam(parIndex.i_EGP_hep)=0;
Optparam(parIndex.i_G0)=11;
Optparam(parIndex.i_G_GTT)=11;

% Select start guess as initial value for the parameters

Optparam(OPTIONS.index_Optpar)=startGuess;

icOrig0=[(11+Optparam(parIndex.i_delta_G_1_11mM))*param(parIndex.i_V_m_liver) ...
    (11+Optparam(parIndex.i_delta_G_1_11mM))*param(parIndex.i_V_m_pancreas)...
   (0+Optparam(parIndex.i_delta_I_1_11mM))*param(parIndex.i_V_m_liver) (0+Optparam(parIndex.i_delta_I_1_11mM))*param(parIndex.i_V_m_pancreas) 0 0 ...
   5.5 0.0000000088];

startCost = costFunction(startGuess);

% Simulated annealing
% Set temperature
OPTIONS.tempstart = 1e1*startCost;                   % InitialTemp
OPTIONS.maxitertemp =   50*length(startGuess);       % Max-iterations per temp
OPTIONS.maxitertemp0 =  200*length(startGuess);      % Max-iterations per temp0

% SBPDsimulate optionsc
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

% Costfunction settings
COSTOPTIONS = [];
COSTOPTIONS.maxnumsteps = 1e8;
COSTOPTIONS.abstol = 1e-10;
COSTOPTIONS.reltol = 1e-10;
COSTOPTIONS.minstep = 0;
COSTOPTIONS.maxstep = inf;

%% Choose data to optimize

%  Cut-off value for chi-2 test

% Calculate number of data points for statistical tests
dataPoints=0;

for i=1:OPTVAR_INDEX(end)
    dataPoints=dataPoints+length(EXPDATA.time{i});
end

CUTOFF = chi2inv(0.95,dataPoints);

file = ['parameterValues' ' ' modelName ' ' datestr(datetime('now')) '.dat'];

FID = fopen(file, 'wt');

N_iter=10;      % Number of iterations for the optimization

for i=1:N_iter
    startGuess=X;
    startCost = costFunction(startGuess);
    disp(['Initial cost with start guess: ' num2str(startCost)]);
    
    [X,FVAL,EXITFLAG] = simannealingSBAOClusteringL(@costFunction,startGuess,OPTIONS);
    
end







