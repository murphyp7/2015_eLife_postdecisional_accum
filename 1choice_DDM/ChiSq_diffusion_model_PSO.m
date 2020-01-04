% Takes input from FIT_diffusion_model_PSO as part of particle swarm
% optimisation - note, multiple parameter sets are passed into this
% function; hence the loop through the fist dimension of params

function [chisq,allRT] = ChiSq_diffusion_model_PSO(params)

global leak
global boundary_decay
global preerror_time

global quants
global O
global n_trials
global maxRT

%% Fixed params
sigma = 0.1;      % noise in decision process
tau = 120;        % evidence decay constant (derived from leaky theta/Pe fit = 120ms)
lambda = 0.00048;  % boundary exponential decay constant (derived from Pe amplitude/latency relationship = 0.00048)
t0_extra = 220;   % time before error commission at which accumulation begins (derived from theta ROC analysis = 220ms)
cutoff = maxRT./1000;     % time beyond which trial is counted as unaware

%% Simulation settings
n_sims = 15000;

for particle = 1:size(params,1);
    %% Free params
    v = params(particle,1);         % drift rate
    eta = params(particle,2);       % drift rate variability
    a = params(particle,3);         % evidence bound
    Ter = params(particle,4);       % non-decision time
    
    %% Beginning trial loop
    RT = [];
    
    for trial = 1:n_sims;
        
        vS = normrnd(v,eta);   % drawing drift rate for current trial
        DV = 0;   % initialising DV to zero
        aS = a;   % intialising bound to a
        
        if strcmp(preerror_time,'OFF')
            t = 0;   % initialising time to zero
        elseif strcmp(preerror_time,'ON')
            t = -t0_extra;
        end
        
        while t/1000<cutoff
            
            t = t+1;   % moving forward one time step
            
            if strcmp(boundary_decay,'ON') && t>0
                aS = a*exp(-lambda*t);   % decreasing the bound in accordance with decay constant lambda
            end
            
            % moving DV by one step (decay will only be enforced if current state of DV is positive)
            sim_rnd = rand(1);
            if sim_rnd < 0.5.*(1+(vS.*(sqrt(0.001)./sigma)));
                if DV(end) > 0 && strcmp(leak,'ON')
                    DV(end+1,1) = DV(end)+(sigma.*sqrt(0.001))-(DV(end)./tau);
                else DV(end+1,1) = DV(end)+(sigma.*sqrt(0.001));
                end
            else
                if DV(end) > 0 && strcmp(leak,'ON')
                    DV(end+1,1) = DV(end)-(sigma.*sqrt(0.001))-(DV(end)./tau);
                else DV(end+1,1) = DV(end)-(sigma.*sqrt(0.001));
                end 
            end
            
            % checking whether bound or cutoff time have been passed
            if DV(end) > aS
                RT(trial) = t./1000+Ter;
                break
            elseif t./1000+Ter >= cutoff;
                RT(trial) = NaN;
                break
            end
            
        end
    end
    
    %% Calculating expected frequencies
    E = [];
    for q = 1:length(quants)
        if q == 1
            E(1,q) = length(find(RT<=quants(1)./1000))./length(RT).*n_trials;
        else E(1,q) = length(find(RT>quants(q-1)./1000 & RT<=quants(q)./1000))./length(RT).*n_trials;
        end
    end
    E(1,end+1) = length(find(RT>quants(end)./1000))./length(RT).*n_trials;
    
    %% Calculating chi-square
    err = [];
    for q = 1:length(O);
        if E(q)>0
            err(q) = ((O(q)-E(q)).^2)./E(q);
        else err(q) = ((O(q)-E(q)).^2)./1;
        end
    end
    chisq(particle,1) = sum(err);
    
    allRT(:,particle) = RT;   % outputting all simulated RT distributions
end


        
