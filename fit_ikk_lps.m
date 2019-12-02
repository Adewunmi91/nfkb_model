function [sse] = fit_ikk_lps(vals, param_idx, varargin)
% param_idx is [n x 2] matrix specifiying full modified parameter index
% vals is altered values
% Fit measured IKK timepoints (continuous stimulation at each of 2 doses


% Fit data columns: 
% { 1 ng LPS } | {100 ng LPS} 
% time |  IKK  | time |  IKK
fit_data = varargin{1};

species_mod = {};
vals = reshape(vals,[numel(vals),1]);

p_mod = [param_idx, vals];
species_mod = {};

scale_factor = 0.036; % max scale of IKK activity (for 100ng/mL dose)


% Model requires TNF in uM -> 1.96e-4 uM = 1ng/mL TNF
options.PULSE_TIME = 45;
options.SIM_TIME = 240;

% LPS input is uM -> 100,000x g/mL (e.g. 1ng/mL = 0.001 uM input)
doses = [1 100];
dose_scale = 1/10000;
names = {'IKK','NFkBn'};

% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, species_mod,options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'LPS',doses(i)*dose_scale},names, p_mod, species_mod,options);
    end
    output = cat(3,output,x);
end

% Define errors from interpolated IKK timecourse
ikk_tp = [0 3 10 25 40 60 90 120];
ikk = squeeze(output(:,strcmp(names,'IKK'),:));
ikk_lo = ikk(ikk_tp+1,1);
ikk_hi = ikk(ikk_tp+1,2);
error_vect = [(ikk_lo - fit_data(:,2)*scale_factor); (ikk_hi - fit_data(:,4)*scale_factor)];
sse = norm(error_vect);

%% Plot simulation results against measured data
if nargin>3
    figure('Position',[500        1084         644         266],'Name','Modeled vs actual LPS IKK response')
    hold on
    plot(ikk,'LineWidth',2)
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(ikk_tp,scale_factor*fit_data(:,[2,4]),'o')
    hold off
    set(gca,'XLim',[0 120])
end
