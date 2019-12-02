for n = [0.3]

sc_sim(1).stim = 'CpG';
sc_sim(1).doses = (10.^(1:.33:3));
sc_sim(1).TLR = [85 1 2e-6];
sc_sim(1).dose_scale = 1/1000;

sc_sim(2).stim = 'LPS';
sc_sim(2).doses =  (10.^(-0.6:.33:1.5));
sc_sim(2).TLR = [35 1 5.25E-05];
sc_sim(2 ).dose_scale = 1/24000;

sc_sim(3).stim = 'Pam3CSK';
sc_sim(3).doses =  (10.^(0:.33:2));
sc_sim(3).TLR = [68 1 1e-6];
sc_sim(3).dose_scale = 1/1500;

% Specify other parameters 
n_sims = 512;
names = {'IKK','NFkBn'};
opt1 = struct;
opt1.DEBUG = 0;
opt1.SIM_TIME = 8*60;

noise_add = n;

for idx = 1:length(sc_sim)
    basal_tlr = sc_sim(idx).TLR(3);
    basal_mydtrif = 0.1;
    all_myd = 10.^(log10(basal_mydtrif) + noise_add*randn(n_sims,1));
    all_trif = 10.^(log10(basal_mydtrif) + noise_add*randn(n_sims,1));
    all_tlr = 10.^(log10(basal_tlr) + noise_add*randn(n_sims,1)); 
    
    % Slice off appropriate structure data
    stim_name = sc_sim(idx).stim;
    doses = sc_sim(idx).doses.*sc_sim(idx).dose_scale;
    tlr_idx = sc_sim(idx).TLR(1:2);
    all_sim = zeros(opt1.SIM_TIME+1, length(names),length(sc_sim(idx).doses),n_sims);
    
    
    
    parfor cell = 1:length(all_myd)
        disp(['Stimulus:', stim_name, '. Cell #', num2str(cell)])
        p_mod = [tlr_idx, all_tlr(cell)];
        s_mod = {'MyD88_off',all_myd(cell); 'TRIF_off',all_trif(cell)};
        % Simulate all doses (only need to equilibrate on first iteration)
        options = opt1;
        output = [];
        for i = 1:length(doses)
            if isempty(output)
                [t,x,simdata] = nfkbSimulate({stim_name,doses(i)},names, p_mod, s_mod,options);
            else
                options.STEADY_STATE = simdata.STEADY_STATE;
                [~,x] = nfkbSimulate({stim_name,doses(i)},names, p_mod, s_mod , options);
            end
            output = cat(3,output,x);
        end
        all_sim(:,:,:,cell) = output;
    end
    sc_sim(idx).data = all_sim;
end

save(['/home/brooks/Data/sc_sim_',num2str(n),'noise.mat'],'sc_sim')
end
    