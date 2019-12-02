function [sse, nfkb, p_mod] = fit_oscillations(vals, param_idx, fit_trajectory, extra_params,varargin)
% param_idx is [n x 3] matrix specifiying full parameter index, and altered value
% vals is altered values

if ~isempty(varargin)
    show_graph = varargin{1};
else
    show_graph = 0;
end

% Fit measured oscillation timepoints in BMDMs
species_mod = {};
vals = reshape(vals,[numel(vals),1]);


% add paired parameters
if ismember(21,param_idx(:,1))
    val_idx = find(param_idx(:,1)==21);
    param_idx = [param_idx; 22 1];
    vals = [vals; vals(val_idx)];
end

if ismember(23,param_idx(:,1))
    val_idx = find(param_idx(:,1)==23);
    param_idx = [param_idx; 24 1];
    vals = [vals; vals(val_idx)];
end

if ismember(25,param_idx(:,1))
    val_idx = find(param_idx(:,1)==25);
    param_idx = [param_idx; 26 1];
    vals = [vals; vals(val_idx)];
end

p_mod = [param_idx, vals];

% Optional: specify additional modified parameters
if nargin>3
    p_mod = [p_mod; extra_params];
end


% Model requires TNF in uM -> 1.96e-4 uM = 1ng/mL TNF
dose_scale = 1.94e-4;
doses = 10.^([1]);
names = {'NFkBn'};
options.SIM_TIME = 400;
% Simulate all doses (only need to equilibrate on first iteration)
output = [];
for i = 1:length(doses)
    if isempty(output)
        [t,x,simdata] = nfkbSimulate({'TNF',doses(i)*dose_scale},names, p_mod, species_mod,options);
    else
        options.STEADY_STATE = simdata.STEADY_STATE;
        [~,x] = nfkbSimulate({'TNF',doses(i)*dose_scale},names, p_mod, species_mod,options);
    end
    output = cat(3,output,x);
end

nfkb = squeeze(output);
nfkb_scaled = (nfkb-min(nfkb))/range(nfkb(:));
fit_trajectory = (fit_trajectory-min(fit_trajectory))/range(fit_trajectory);

% Define errors from a "simple" trajectory
error_vect = nfkb_scaled(1:length(fit_trajectory)) - fit_trajectory';
sse = norm(error_vect);


if show_graph
    figure('Position',[500        1084         644         266],'Name','Modeled vs actual TNF response')
    hold on
    plot(nfkb_scaled,'LineWidth',2,'Color',[0.2902 0.2980 0.3098])
    plot(fit_trajectory,'Color',[0.4627 0.7059 0.7961],'LineWidth',2)
    hold off
    set(gca,'XLim',[0 400],'YLim',[-0.01 1.05])
end
