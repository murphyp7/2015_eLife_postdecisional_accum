% function [best_params,chisq] = FIT_diffusion_model_PSO_SS(subj)
%
% Fits the one-choice diffusion model to single-subject error signalling
% RT distributions, using particle swarm optimzation.
%
% Input: 'subj' = string variable representing participant number
%
% Outputs: 'base_params' = vector of best-fitting diffusion parameters
%          'chisq' = chi-squared value minimized by PSO procedure
%          'tr' = chi-sq minimum at every iteration
%          'te' = number of iterations taken to optimize
%
% Peter Murphy. UKE Hamburg, 04/01/2020

function [best_params,chisq,tr,te] = FIT_diffusion_model_PSO_SS(subj)

%% Assigning global variables to pass to model fitting
global leak
global boundary_decay
global preerror_time

global quants
global O
global n_trials
global maxRT

%% Adding required paths
addpath(genpath('/home/murphy/EAT_DDM'))
addpath(genpath('/home/murphy/matlab/Imported_toolboxes'))

%% Define saving directory
savepath = '/home/murphy/EAT_DDM/OUTPUT/';

%% Choosing whether to impose EEG-informed model constraints
leak = 'OFF';              % accumulator leak
boundary_decay = 'OFF';    % collapsing bound
preerror_time = 'OFF';     % time prior to error commission at which accumulation onsets

%% Loading behavioural data
subj = num2str(subj);
loadpath = '/home/murphy/EAT_DDM/RT_dists/';
load([loadpath,subj,'_AwareRT_dist.mat'])

%% Seed random number generator
seed = round(sum(100*clock)); %never the same seed
rand('state', seed);

%% Calculating essential measures for fitting
n_trials = length(Aware_RTs)+n_unawares;   % total number of error trials for this participant
awareness = length(Aware_RTs)./n_trials;   % percentage of those trials on which awareness was achieved
maxRT = max(Aware_RTs);

% Calculating observed RT quantiles and trial frequencies therein
quants = quantile(Aware_RTs,[0.1 0.3 0.5 0.7 0.9]);   % getting error signalling RT quantiles
O = [awareness.*n_trials.*0.1 awareness.*n_trials.*0.2 awareness.*n_trials.*0.2 awareness.*n_trials.*0.2 awareness.*n_trials.*0.2 awareness.*n_trials.*0.1];   % getting observed frequencies for each quantile

%% Setting PSO parameters
% Defining maximum particle velocity for each parameter
mv = [0.1 0.05 0.03 0.03];

% Defining parameter ranges
par_range = [0.01 2;...  % v
            0.001 1.2;... % eta
            0.03 0.8;... % a
            0.01 0.6];   % Ter
        
% Defining PSO options (see pso_Trelea_vectorized.m for details)
P(1)=0;  P(2)=2000;  P(3)=50;  P(4:13)=[1.6 1.9 0.9 0.4 1500 0.01 250 NaN 0 1];

% Defining initial seed values (1 row per particle)
PSOseedValue = [...
    0.8 0.3 0.13 0.18; ...
    1 0.5 0.18 0.13; ...
    0.4 0.3 0.09 0.19; ...
    0.74 0.12 0.25 0.1; ...
    0.15 0.1 0.1 0.1; ...
    1.1 0.11 0.45 0.3; ...
    0.9 0.65 0.3 0.34; ...
    0.3 0.1 0.11 0.3; ...
    0.4 0.4 0.5 0.2; ...
    1.3 0.8 0.22 0.24; ...
    0.36 0.19 0.27 0.26; ...
    0.4 0.005 0.14 0.25; ...
    0.89 0.06 0.17 0.1; ...
    0.49 0.10 0.10 0.07; ...
    0.51 0.3 0.13 0.37; ...
    0.6 0.21 0.08 0.4; ...
    0.7 0.1 0.34 0.24; ...
    0.91 0.64 0.39 0.22; ...
    0.075 0.2 0.05 0.24; ...
    0.8 0.3 0.21 0.21; ...
    0.66 0.06 0.22 0.17; ...
    1.05 0.25 0.43 0.19; ...
    0.65 0.25 0.13 0.31; ...
    0.3 0.01 0.12 0.05; ...
    0.85 0.4 0.14 0.45];

PSOseedValue(size(PSOseedValue,1)+1:P(3),1:size(PSOseedValue,2)) = normmat(rand([5,size(PSOseedValue,2)]),par_range',1);  % including an extra 5 particles with randomized initial positions

%% Running PSO routine
[output,tr,te,simRT] = pso_Trelea_vectorized('ChiSq_diffusion_model_PSO',4,mv,par_range,0,P,'goplotpso',PSOseedValue);

%% Saving output
best_params = output(1:end-1);
chisq = output(end);

save([savepath,subj,'_DDM_output_EEG_',leak,'.mat'],'best_params','chisq','tr','te','n_trials','awareness','maxRT','quants','O','Aware_RTs','simRT')
