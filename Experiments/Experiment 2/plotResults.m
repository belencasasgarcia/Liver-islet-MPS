% Plot model simulations with uncertainty based on the estimated parameter
% values

clear
clc
close all
format long
format compact

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
global statesIndex
global statesNames

%% Load design parameters for the 2-OC system

MPSDesignParameters

%% Load and plot experimental data

loadExperimentalData

% Define number of media exchanges (time intervals of 48 h) to simulate

N_exchanges=7;

modelName ='MPSmodel';

optModel = IQMmodel(strcat(modelName,'.txt')); 
IQMmakeMEXmodel(optModel,modelName); 
[pNames, param] = IQMparameters(optModel);
icOrig0 = IQMinitialconditions(optModel);
[statesNames,ODEs,initialConditions] = IQMstates(optModel);

SIMTIME=[0:0.01:336];

% Load parameter indexes
loadParameterIndexes

% Load state indexes
loadStateIndexes

% Load optimal parameter values

file = ['parameterValuesExperiment2.dat'];

tmpHold = load(file);

if isempty(tmpHold)
    disp('Empty parameter dataset');
else
    X = getSpreadedParam(tmpHold(:,2:end));
    X_opt = tmpHold(end,2:end);
    COST_OPT=min(tmpHold(:,1));
    plotFunctionUncertainty(X,modelName,parIndex,X_opt,N_exchanges,statesIndex)
end




