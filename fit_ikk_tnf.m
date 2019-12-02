function [sse,ikk,full_params] = fit_ikk_tnf(vals, param_idx, varargin)
% param_idx is [n x 2] matrix specifiying full modified parameter index
% vals is altered values
% Fit measured IKK timepoints (45 min pulse, population-level experiment)
fit = varargin{1};

species_mod = {};
vals = reshape(vals,[numel(vals),1]);


% add paired parameters (TNF degradation)
if ismember(58,param_idx(:,1))
    val_idx = find(param_idx(:,1)==58);
    param_idx = [param_idx; 61 1; 64 1];
    vals = [vals; vals(val_idx); vals(val_idx)];
end
% TNFR snythesis/degradation
if ismember(54,param_idx(:,1))
    val_idx = find(param_idx(:,1)==54);
    param_idx = [param_idx; 55 1];
    vals = [vals;  vals(val_idx)/(3.45e-4)];
end

full_params = [param_idx, vals];

scale_factor = 0.036; % max scale of IKK activity - Werner et al (2008) put peak at roughly 0.035-0.04

% Model requires TNF in uM -> 1.96e-4 uM = 1ng/mL TNF
dose_scale = 1/52000; % TNF molecular weight (trimer) is approx 52kDa
dose = 3.16;
names = {'TNF','TNFR','TNFR_TNF','C1','IKK','NFkBn'};
opt1.PULSE_TIME = 45;
opt1.SIM_TIME = 180;
[t,x] = nfkbSimulate({'TNF',dose*dose_scale},names, full_params, species_mod, opt1);

% Define errors from interpolated IKK timecourse
ikk = x(:,strcmp(names,'IKK'));
fit_scaled = fit.ikk/max(fit.ikk)*scale_factor;
% Construct error vector: allow some "slop" in timecourse collection
slop = 1;
error_vect = zeros(size(fit_scaled));
for i = 1:length(error_vect)
    t_range = fit.t(i)+1+(-slop:slop);
    t_range(t_range<1) = [];
    error_vect(i) = min(abs(ikk(t_range) - fit_scaled(i)));
end
sse = norm(error_vect);


if nargin>3
    figure('Position',[500        1084         644         266],'Name','Modeled vs actual TNF IKK response')
    hold on
    plot(t,ikk,'LineWidth',2,'Color',[0.2902 0.2980 0.3098])
    plot(fit.t,fit_scaled,':o','Color',[0.4627 0.7059 0.7961],'LineWidth',2)
    hold off
    set(gca,'XLim',[0 180])
end
