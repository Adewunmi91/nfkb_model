% - - - - - - - CpG SIMULATION (match doses)- - - - - - - -
p_mod = [
%85 1 0.5e-6
%86 1 1e-4
%87 1 6e-4
%88 1 4e-4 % Degradation of TLR9_N (sets excess amt based on ratio w/ C-terminius degrdataion, 4e-4)
%92 1 3 % Degradation rate via TLR9_N
93 1 1.6e-3 % Constitutive degradation (unbound = 4e-4)
89 1 0.028
];

doses = [10 33 66 100 200 330 1000];
dose_scale = 1/1000; % Convert to uM (from nM)
names = {'CpG','CpG_en','TLR9','TLR9_CpG','TLR9_N','MyD88','TRAF6','IKK','NFkBn'};
options = struct;
options.DEBUG = 1;
options.SIM_TIME = 12*60;
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

plotSpecies(output,names,doses);