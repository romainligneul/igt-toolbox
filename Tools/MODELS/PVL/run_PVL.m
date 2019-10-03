function [A R] = run_EV(A, data)
%%%% PVL model / Evolution function
% This model is well described in Steingroever, Wetzels, Wagenmakers 2016
% 'Bayes Factors for Reinforcement-Learning Models of the IowaGambling Task'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% import options from launcher
options = A.fit.options;

%% define functions used by the model and parameter transformation
% evolution function
evof = @e_PVL;
% observation function
obsf = @o_PVL;

%% define dimensions of the model and priors over parameters
%
dim.n = 4;
dim.n_phi = 1;
dim.n_theta = 3;

switch A.fit.priors.type

    case 'flat' % use strictly bounded priors (flatness imperfect)

        priors.muPhi = zeros(1,dim.n_phi);
        priors.SigmaPhi = 3e0*eye(dim.n_phi);
        priors.muTheta = zeros(1,dim.n_theta);
        priors.SigmaTheta = 3e0*eye(dim.n_theta);
        priors.muX0(1:4,1) = [0 0 0 0];


        Traw = @(x) x;
        Tsig = @sig;
        Texp = @exp;
        Tsig0to5 = @(x) sig(x)*5;
        TsigMin2to2 = @(x) sig(x)*4-2;

    case 'informed' % priors from first pass

        priors.muPhi = -2.1961;
        priors.SigmaPhi = 0.7915^2;
        priors.muTheta = [0.9560    1.6392    1.4903];
        priors.SigmaTheta = [0.955983560352368,0,0;0,1.63919054524668,0;0,0,1.49026435762038;].^2;
        priors.muX0(1:4,1) = [0 0 0 0];

        Traw = @(x) x;
        Tsig = @sig;
        Texp = @exp;
        Tsig0to5 = @(x) sig(x)*5;
        TsigMin2to2 = @(x) sig(x)*4-2;

    case 'shrinkage'

        priors = []; % default settings of the VBA toolbox
        priors.muX0(1:4,1) = [0 0 0 0];

        Traw = @(x) x;
        Tsig = Traw;
        Texp = Traw;
        Tsig0to5 = Traw;
        Tsig0to2 = Traw;
        TsigMin2to2 = Traw;

end

%% transformation of parameters
% inF = evolution
options.inF.param_name = {'Sensitivity (value)', 'Loss aversion (value)', 'Inverse decay (value)'};
options.inF.param_transform = {Tsig Tsig0to5 Tsig};
% inF = evolution (learning)
options.inG.param_name = {'consistency'};
options.inG.param_transform = {Tsig0to5};

%% build hidden states of interest.
%%%% this is a (boring vector) defining the initial values of the hidden
%%%% states inferred by the experimenter.
hs.decks_ind = [1 2 3 4]; % <=> A,B,C,D will be in hidden states x(1), x(2), etc.
hs.initval = 0;% = ones(size(hs.decks_ind))*0; % hidden states are initialized at 0.25
% assign this additional information to optional inputs
options.inF.hs = hs;
options.inG.hs = hs;
% skip observations?
options.skipf = zeros(1,100); % apply identity mapping from x0 to x1.
options.skipf(1) = 1;         % first evolution not taken into account (null trial)

A.fit.options = options;
A.fit.evof = evof;
A.fit.obsf = obsf;

%% perform subject-wise or hierarchical fitting

for s = 1:length(A.fit.u)

    disp(['%%%%%%%%%% SUBJECT ' num2str(s) ' %%%%%%%%%%%']);

    % uses subject-specific priors if available (method second-pass)
    try
        options.priors = priors{s};
    catch
        options.priors = priors;
    end

    % perform fit
    if isfield(A.fit, 'fminunc') &&  A.fit.fminunc == 1
        [posterior{s}, out{s}] = VBA_fminunc_wrapper(A.fit.y{s},A.fit.u{s},evof,obsf,dim,options);
    else
        [posterior{s}, out{s}] = VBA_NLStateSpaceModel(A.fit.y{s},A.fit.u{s},evof,obsf,dim,options);
    end

    % log info
    R.GoF(s,1) =  out{s}.F;
    R.GoF(s,2) =  out{s}.fit.BIC;
    R.GoF(s,3) =  out{s}.fit.AIC;
    R.GoF(s,4) =  out{s}.fit.LL;

    % evolution parameters
    for pp = 1:length(posterior{s}.muTheta)
        R.theta(s,pp) = options.inF.param_transform{pp}(posterior{s}.muTheta(pp));
    end;
    R.rawMuTheta(s,:) = posterior{s}.muTheta;
    R.rawSigmaTheta(s,:,:) = posterior{s}.SigmaTheta;

    % observation parameters
    for pp = 1:length(posterior{s}.muPhi)
        R.phi(s,pp) = options.inG.param_transform{pp}(posterior{s}.muPhi(pp));
    end;
    R.rawMuPhi(s,:) = posterior{s}.muPhi;
    R.rawSigmaPhi(s,:,:) = posterior{s}.SigmaPhi;

    % model dynamics
    R.hidden_states(s,:,:) = posterior{s}.muX;
    for t = 1:size(posterior{s}.muX,2)
        R.actual_choices(s,t) = find(A.fit.y{s}(:,t));
        try
            R.predicted_choices(s,t) = find(out{s}.suffStat.gx(:,t)==max(out{s}.suffStat.gx(:,t)));
        catch
            R.predicted_choices(s,t) = NaN;
        end
    end
    R.choice_prob(s,:,:) = out{s}.suffStat.gx;

    % stop analysis at maxsujects if required
    if s == A.fit.maxsubjects
        break
    end

end


%% simulate & recover

if A.simulate_and_recover && A.fit.fminunc==0

    %%% prepare feedbacks for simulation analysis
    for d = 1:4
        deck_fb{d} = [];
    end
    for s = 1:length(data)
        for d = 1:4
            deck_fb{d} = [deck_fb{d};[data{s}.win(data{s}.deck==d) data{s}.lose(data{s}.deck==d)]/A.fit.divide_feeback];
        end
    end
    for d = 1:4
        ranges(d,:) = [min(deck_fb{d}(:,1)) max(deck_fb{d}(:,1)) min(deck_fb{d}(deck_fb{d}(:,2)>0,2)) max(deck_fb{d}(:,2)) ];
    end

    fb.h_fname = @feedback_IOWA;
    fb.inH = deck_fb;
    % allocate input for feedback simulation
    fb.indfb = [3 4];
    fb.indy = [7:10];

    for s=1:length(A.fit.u)

        % set-up simulation
        theta = posterior{s}.muTheta;
        phi = posterior{s}.muPhi;
        u = out{s}.u;
        u(7:10,:) = 0;
        simopts = out{s}.options;
        simopts.DisplayWin = 0;
        x0 = posterior{s}.muX0;
        n_t = 100;

        % recover parameters (sometimes the simulation fails due to
        % aberrant series of outcomes which exceed the numerical capacities
        % of matlab, ie when x>exp(710), then try until ok)
        err_count = 0; count = 0;
        while err_count == count
            try
                [sim_y,x,~,~,~,sim_u] = simulateNLSS(n_t,evof,obsf,theta,phi,u,Inf,Inf,simopts,x0,fb);
                [sim_post sim_out] = VBA_NLStateSpaceModel(sim_y,sim_u,evof,obsf,simopts.dim,simopts);
            catch
                err_count = err_count + 1;
            end
            count = count+1;
        end

        % log info of the re-fit
        R.simulation.theta_recovered(s,:) = sim_post.muTheta;
        R.simulation.phi_recovered(s,:) = sim_post.muPhi;
        for t = 1:size(posterior{s}.muX,2)
            R.simulation.simulated_choices(s,t) = find(sim_y(:,t));
            try
                R.simulation.predicted_choices(s,t) = find(out{s}.suffStat.gx(:,t)==max(out{s}.suffStat.gx(:,t)));
            catch
                R.simulation.predicted_choices(s,t) = NaN;
            end
        end
        R.simulation.simulated_gains(s,:) = sim_u(3,:);
        R.simulation.simulated_loss(s,:) = sim_u(4,:);
        R.simulation.choice_prob(s,:,:) = sim_out.suffStat.gx;

        % log GoF
        R.simulation.GoF(s,1) =  sim_out.F;
        R.simulation.GoF(s,2) =  sim_out.fit.BIC;
        R.simulation.GoF(s,3) =  sim_out.fit.AIC;
        R.simulation.GoF(s,4) =  sim_out.fit.LL;

        % stop analysis at maxsujects if required
        if ~isempty(A.fit.maxsubjects) && s == A.fit.maxsubjects
            break
        end

    end
end

% decide numbering
dt = 1;
while exist([A.output.dir char(obsf) '_' char(evof) '_' date '_' num2str(dt)])
    dt = dt+1;
end;
if A.fit.fminunc == 1
    R.output_subdir=[A.output.dir 'fminunc_' char(obsf) '_' char(evof) '_' date '_' num2str(dt)];
else
    R.output_subdir=[A.output.dir char(obsf) '_' char(evof) '_' date '_' num2str(dt)];
end

% make and fill dir
mkdir(R.output_subdir);
if A.complete_save==1
    save([R.output_subdir '/fitted_model'], 'out', 'posterior', 'R', 'A');
elseif A.complete_save==-1
    save([R.output_subdir '/fitted_model'], 'R');
else
    save([R.output_subdir '/fitted_model'], 'R', 'A');
end
% save script in dir
scriptname  = mfilename('fullpath');
copyfile([scriptname '.m'],[R.output_subdir '/generating_script.m'])

end
