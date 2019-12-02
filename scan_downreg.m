 %% Script to explore frequency/steady state space enforced by different transport rates.

% Define multipliers/discret delays/halflives of appropriate species
param.stim_lvl = 10.^linspace(-3, -1, 20);
param.stim_time = linspace(60, 240, 10);
param.downreg = log(2)./linspace(10,100,10);
[X,Y,Z] = ndgrid(param.stim_lvl,param.stim_time,param.downreg);
all_param = [X(:),Y(:),Z(:)];

% Default runout is 24 hrs, or 1441 min. Downsample slightly (2 min timepoints - 721 timepoints)
all_nfkb_n = single(zeros(721,size(all_param,1)));
all_ikk = single(zeros(721,size(all_param,1)));



parfor idx = 1:size(all_param,1)
    disp(['parameter alterations: [', num2str(all_param(idx,:)),']'])  
    p_mod = [3 1 all_param(idx,3)
             2 1 0]; 
    names = {'IKK','NFkBn'};
    options = struct;
    options.SIM_TIME = 24*60;
    options.PULSE_TIME = all_param(idx,2);

    % Simulate all doses (only need to equilibrate on first iteration)
    output = [];  
    [~,sim_x] = nfkbSimulate({'stim',all_param(idx,1)},names, p_mod, {},options);
    
    % Record output - NFkB trajectories, and features (min steady state, max steady state,
    all_nfkb_n(:,idx) = single(squeeze(sim_x(1:2:end,2))');
    all_ikk(:,idx) = single(squeeze(sim_x(1:2:end,1))');
end


% Reshape things for easier indexing 
all_nfkb_n = reshape(all_nfkb_n,[size(all_nfkb_n,1),length(param.stim_lvl),length(param.stim_time),length(param.downreg)]);
all_ikk = reshape(all_ikk,[size(all_ikk,1),length(param.stim_lvl),length(param.stim_time),length(param.downreg)]);

save('/home/brooks/Data/scan_downreg.mat','all_nfkb','all_ikk','param')



%% 2nd form of downregulation: sequestration (represented here as degradadation of NfkB)
% Define multipliers/discret delays/halflives of appropriate species
param.stim_lvl = 10.^linspace(-1, 1, 20);
param.stim_time = linspace(60, 240, 10);
param.downreg = log(2)./linspace(1,10,10);
 [X,Y,Z] = ndgrid(param.stim_lvl,param.stim_time,param.downreg);
all_param = [X(:),Y(:),Z(:)];

% Default runout is 24 hrs, or 1441 min. Downsample slightly (2 min timepoints - 721 timepoints)
all_nfkb_n = single(zeros(721,size(all_param,1)));
all_nfkb_c = single(zeros(721,size(all_param,1)));

all_ikk = single(zeros(721,size(all_param,1)));


parfor idx = 1:size(all_param,1)
    
    disp(['parameter alterations: [', num2str(all_param(idx,:)),']'])  
    names = {'IKK','NFkBn','IkBaNFkBn','NFkB','IkBaNFkB','IKKIkBaNFkB'};
    options = struct;

    % Simulate 1st period w/o degradataion
    p_mod = [];

    options.SIM_TIME = all_param(idx,2);
    [~,x1,~,end_state] = nfkbSimulate({'stim',all_param(idx,1)},names, p_mod, {},options);

    % Pass final values as "steady state" (i.e. initial values) and simulate next 23 hrs
    p_mod = [95 1 all_param(idx,3)];
    options.SIM_TIME = 24*60 - all_param(idx,2) -1;
    options.STEADY_STATE = end_state;
    [~,x2] = nfkbSimulate({'stim',all_param(idx,1)},names, p_mod, {},options);

    sim_x = [x1;x2];
    
    % Record output - NFkB trajectories, and features (min steady state, max steady state,
    all_nfkb_n(:,idx) = single(squeeze(sim_x(1:2:end,2))' + squeeze(sim_x(1:2:end,3))');
    all_nfkb_c(:,idx) = single(squeeze(sim_x(1:2:end,4))' + squeeze(sim_x(1:2:end,6))' + squeeze(sim_x(1:2:end,5))');

    all_ikk(:,idx) = single(squeeze(sim_x(1:2:end,1))');

end

% Reshape things for easier indexing 
all_nfkb_n = reshape(all_nfkb_n,[size(all_nfkb_n,1),length(param.stim_lvl),length(param.stim_time),length(param.downreg)]);
all_nfkb_c = reshape(all_nfkb_c,[size(all_nfkb_c,1),length(param.stim_lvl),length(param.stim_time),length(param.downreg)]);
all_ikk = reshape(all_ikk,[size(all_ikk,1),length(param.stim_lvl),length(param.stim_time),length(param.downreg)]);

save('/home/brooks/Data/scan_downreg_sequester.mat','all_nfkb_c','all_nfkb_n','all_ikk','param')