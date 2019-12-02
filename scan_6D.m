e%% Script to explore frequency/steady state space enforced by different transport rates.

% Define multipliers/discret delays/halflives of appropriate species
param.nfkb = [1.8 0.6 0.2];
param.ikba = [0.27 0.09 0.3];
param.dd1 = [11 14 17 20 23];
param.dd2 = [1 4 7];
param.ikba_hl = [4 6 8];
param.ikbat_hl = [12 16 20 24 28];
[X,Y,Z,P,D,Q,M] = ndgrid(param.nfkb,param.ikba,param.dd1,param.dd2,param.ikba_hl,param.ikbat_hl);
all_param = [X(:),Y(:),Z(:),P(:),D(:),Q(:),M(:)];
doses = 10.^(linspace(-5,4,48));

% Default runout is 72 hrs, or 4321 min. Downsample slightly (2 min timepoints - 2161 timepoints)
all_nfkb = single(zeros(length(doses),2161,size(all_param,1)));
all_ss = nan(length(doses),2,size(all_param,1));
all_freq = nan(length(doses),1,size(all_param,1));

%%
parfor idx = 1:size(all_param,1)
    disp(['parameter alterations: [', num2str(all_param(idx,:)),']'])  
    p_mod = [
            % IkBa transport: rxns #9 and #11
            9   1   all_param(idx,2)
            11  1   all_param(idx,2)/2*3.5
            % NFkB transport: rxns #10 and #12
            10  1   all_param(idx,1)
            12  1   all_param(idx,1)/50*3.5
            % discrete delay 1 (mRNA)
            6   4   all_param(idx,3)
            % discrete delay 2 (protein)
            8   2   all_param(idx,4)
            % IkBa halflife
            15  1   log(2)/all_param(idx,5)
            16  1   log(2)/all_param(idx,5)
            % IkBat halflife
            7   1   log(2)/all_param(idx,6)
            ];

        % Run dose response scan - no interpolation, no graphing.
        options = struct;
        [simdata, features] = stimscan(doses,'ModParams',p_mod,'Plots',0,'Parallel',0);

        % Record output - NFkB trajectories, and features (min steady state, max steady state,
        all_nfkb(:,:,idx) = single(squeeze(simdata(1:2:end,2,:))');
        all_ss(:,:,idx) = features.steady_state';
        all_freq(:,:,idx) = features.periods';
end

%%
% Reshape things for easier indexing 
all_nfkb = reshape(all_nfkb,[size(all_nfkb,1),size(all_nfkb,2),length(param.nfkb),length(param.ikba),length(param.dd1),length(param.dd2),length(param.ikba_hl),length(param.ikbat_hl)]);
all_ss = reshape(all_ss,[size(all_ss,1),size(all_ss,2),length(param.nfkb),length(param.ikba),length(param.dd1),length(param.dd2),length(param.ikba_hl),length(param.ikbat_hl)]);
all_freq = reshape(all_freq,[size(all_freq,1),size(all_freq,2),length(param.nfkb),length(param.ikba),length(param.dd1),length(param.dd2),length(param.ikba_hl),length(param.ikbat_hl)]);
save('/home/brooks/Data/scan_7D.mat','all_nfkb','all_ss','all_freq','all_param','doses','param')

