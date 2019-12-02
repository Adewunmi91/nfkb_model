%% script to explore frequency/steady state space enforced by different transport rates.
mult_ikba = 10.^linspace(-1,1,13);
mult_nfkb = 10.^linspace(-1,1,13);
mult_complex = 10.^linspace(-1,1,13);
doses = 10.^(linspace(-5,4,48));

% Make meshgrid to get all screen combinations (rearrange slightly to account for [r,c] ordering)
[X,Y,Z] = meshgrid(mult_nfkb,mult_ikba,mult_complex);
all_mult = [X(:),Y(:),Z(:)];
 
% Default runout is 72 hrs, or 4321 min. Downsample slightly (2 min timepoints - 2161 timepoints)
all_nfkb = single(zeros(length(doses),2161,size(all_mult,1)));
all_ss = nan(length(doses),2,size(all_mult,1));
all_freq = nan(length(doses),1,size(all_mult,1));

%%
parfor idx = 1:size(all_mult,1)
    disp(['Current multipliers: [',num2str(all_mult(idx,2)),'x, ',num2str(all_mult(idx,1)),'x, ',num2str(all_mult(idx,3)),'x ]'])
            p_mod = [
                % IkBa transport: rxns #9 and #11
                9   1   0.09*all_mult(idx,2)
                11  1   0.09/2*3.5*all_mult(idx,2)
                % NFkB transport: rxns #10 and #12
                10  1   0.6*all_mult(idx,1)
                12  1   0.6/50*3.5*all_mult(idx,1)
                % IkBa-NFkB transport: rxn #14
                14  1   0.828*all_mult(idx,3)
                ];
            
            % Run dose response scan - no interpolation, no graphing.
            options = struct;
            [simdata, features] = stimscan_series(doses, p_mod, {},options, 0,0);
            
            % Record output - NFkB trajectories, and features (min steady state, max steady state,
            all_nfkb(:,:,idx) = single(squeeze(simdata(1:2:end,2,:))');
            all_ss(:,:,idx) = features.steady_state';
            all_freq(:,:,idx) = features.periods';
end

%%
% Reshape things for easier indexing 
all_nfkb = reshape(all_nfkb,[size(all_nfkb,1),size(all_nfkb,2),length(mult_ikba),length(mult_nfkb),length(mult_complex)]);
all_ss = reshape(all_ss,[size(all_ss,1),size(all_ss,2),length(mult_ikba),length(mult_nfkb),length(mult_complex)]);
all_freq = reshape(all_freq,[size(all_freq,1),size(all_freq,2),length(mult_ikba),length(mult_nfkb),length(mult_complex)]);
save('/home/brooks/Data/transport_scan.mat','all_nfkb','all_ss','all_freq','all_mult')

