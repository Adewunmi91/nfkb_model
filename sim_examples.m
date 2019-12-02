% TNF SIMULATION
%p_mod = [idx, all_x(:,idx(1))];
dose_scale = 1/5200; % Dose scaling from Shannon et al (2007)
doses = [0.3 3.3 33];
names = {'TNF','TNFR','TNFR_TNF','C1','IkBat','IkBa','IkBan','IKK','NFkBn','polyIC', 'Pam3CSK'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;
% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'TNF',doses(i)*dose_scale},names, [], {},options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'TNF',doses(i)*dose_scale},names, [], {},options);
    end
    output = cat(3,output,x);
end

ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));


%% - - - - - - - LPS SIMULATION (match doses)- - - - - - - -
p_mod = [];
doses = [0.33 3.3 33];
dose_scale = 1/24000; % LPS molecular weight estimated between 10KDa and 100KDa
names = {'TLR4','TLR4LPS','TLR4LPSen','TRIF','MyD88','TRAF6','IKK','IkBat','NFkBn'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 48*60;
% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, {},options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, {},options);
    end
    output = cat(3,output,x);
end

%plotSpecies(output,names,doses);

ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));

%% - - - - - - - CpG SIMULATION (match doses)- - - - - - - -
p_mod = [
%88 1 4e-4 % Degradation of TLR9_N (sets excess amt based on ratio w/ C-terminius degrdataion, 4e-4)
%92 1 3 % Degradation rate via TLR9_N
%93 1 1.6e-3 % Constitutive degradation (unbound = 4e-4)
%89 1 0.028
];
doses = [10 33 100 330];
dose_scale = 1/1000; % Convert to uM (from nM)
names = {'CpG','CpG_en','TLR9','TLR9_CpG','TLR9_N','MyD88','TRAF6','IKK','NFkBn'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;
% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'CpG',doses(i)*dose_scale},names, p_mod, {},options);
    end
    output = cat(3,output,x);
end


ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));

%% - - - - - - - poly(I:C) SIMULATION (match doses)- - - - - - - -
p_mod = [];
doses = 1000*[3.3 10 33 100];
dose_scale = 1/5e6; % Convert to uM. PolyI:C molecular weight: 1000KDa(+)
names = {'polyIC','polyIC_en','TLR3','TLR3_polyIC','TRIF','TRAF6','IKK','NFkBn'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;
% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'polyIC',doses(i)*dose_scale},names, p_mod, {},options);
    end
    output = cat(3,output,x);
end

ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));

%% - - - - - - - Pam3CSK SIMULATION (match doses)- - - - - - - -
p_mod = [];
doses = [1 3.3 10 33];
dose_scale = 1/1500; % Convert to uM. Pam3CSK molecular weight: 1.5KDa
names = {'CD14',  'Pam3CSK', 'CD14_P3CSK','TLR2','TLR2_P3CSK','MyD88','TRAF6','TAK1','IKK','NFkBn'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 8*60;

% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'Pam3CSK',doses(i)*dose_scale},names, p_mod, {},options);
    end
    output = cat(3,output,x);
end

ikk_curves = squeeze(output(:,strcmp(names,'IKK'),:));
nfkb_curves = squeeze(output(:,strcmp(names,'NFkBn'),:));

