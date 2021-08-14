% Plot model simulations with uncertainty based on the estimated parameter
% values

clear
clc
close all
format long
format compact

adjustSigma=boolean(1);  %Parameter to adjust insulin secretion capacity of beta cells
                         % adjustSigma=0: Initial insulin secretion capacity in the
                         % MPS (Fig. S3)
                         % adjustSigma=1: Corrected insulin secretion capacity (Fig. 8)
                
% Declare global variables

global EXPDATA_G
global EXPDATA_I
global EXPDATA_G_hep
global modelName
global icOrig0
global pNames
global param
global OPTIONS
global COSTOPTIONS
global EXPDATA
global SIMTIME      % Time axis for simulation
global OPTIME       % Time axis for optimization
global N_VAR        % Number of variables to optimize
global VAR_INDEX    % Indexes of variables to be optimized
global CUTOFF
global dataPoints
global FID
global Optparam
global parIndex

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

plotExperimentalData(EXPDATA)

modelName ='MPSModel';

optModel = IQMmodel(strcat(modelName,'.txt')); 
IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);

SIMTIME=[0:0.01:384];

% Load parameter indexes

loadParameterIndexes

% SBPDsimulate options
OPTIONS.outputFunction ='';
OPTIONS.silent = 0;
OPTIONS.tempend =       0.1;                     % EndTemp
OPTIONS.tempfactor =    0.1;                     % tempfactor
OPTIONS.maxtime =       120;                     % Max-time
OPTIONS.tolx =          1e-10;                   % TolX
OPTIONS.tolfun =        1e-10;                   % Tolfun
OPTIONS.MaxRestartPoints =  0;                   % Number of parallel valleys which are searched through

% Extrapolate model parameters to translate to humans
 
V_hep_scaled = V_hep*1e5;
V_islets_scaled = V_islets*1e5;
V_m_blood=0.58*5.1;
V_m_hep_scaled = V_m_blood/2;
V_m_islets_scaled = V_m_blood/2;
Total_V_m=V_m_hep_scaled+V_m_islets_scaled;

%Recirculation time in the human body is about 5 minutes=5/60 hours
r_time=5/60;
Q_scaled=(Total_V_m)/r_time;

% Optimal parameter values
file = ['parameterValuesExperiment1Gttd13.dat'];


tmpHold = load(file);


if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    
    parMPS = getSpreadedParam(tmpHold(:,2:end));
    parMPS_opt = tmpHold(end,2:end);
    
    % Extrapolate parameters for human translation
    
    [parHuman V_islets_Human]=extrapolateToHuman(parMPS,V_hep,V_islets,parIndex,adjustSigma);
   
    plotFunctionUncertainty(parHuman,modelName,parIndex,parMPS_opt,V_islets_Human)

end




