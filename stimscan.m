function [sim_output, features, figs, regimes] = stimscan(doses, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [sim_output, features, figs] = stimscan(doses, p_mod,species_mod, options, show_graphs, do_interp)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% STIMSCAN performs a dose response curve @ the level of IKK activation.
% 
% INPUTS (required):
% doses    range of doses (of 'stim' or 'TAK1') to simulate
% 
% INPUTS (optional) - passed as name-argument pairs
% 'ModParams'     change model parameters - [n x 3] matrix w/ row specifying: [idx1 , idx2 , val] (e.g. [6  2   0.003])
%                     (defaults to empty, or no modification)
% 'ModSpecies'    change initial amount of species (e.g. NFkB) - [n x 2] cell matrix w/ rows of species name, then value
%                     (defaults to empty, or no modification)
% 'Options'       change simulation options (structure; see nfkbSimulate for fields)
% 'Parallel'      boolean flag sets whether to simulate each dose in parallel, or in series. (default: parallel)
% 'Plots'         boolean flag to show plots (will also cause output feature to be interpolated)
% 'Stimulus'      set input to generalized non-cooperative stimulus ('stim'), or cooperative TAK1. (Default: 'stim')
% 'SimData'       if provided, skip to feature-finding and graphing
%
% OUTPUTS:
% sim_output     NFkB and IKK curves (sim_output)
% features       oscillatory features (period and steady state)
% figs           figure handles
% regimes        boundaries for "biological cutoffs" (e.g stable oscillation, damped oscillation, etc)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: doses.
addRequired(p,'doses',@isnumeric);

% Optional: parameter modifier matrix (each row specifies [idx1 idx2 val])
valid_pmod = @(x) assert((isempty(x)||(size(x,2)==3))&&(isnumeric(x)), ...
    'Specify parameter modification using must be of size [n x 3]');
addParameter(p,'ModParams',[],valid_pmod);
% Optional: species modifier cell (changes species steady state - name/value pairs
valid_species = @(x) assert((size(x,2)==2)&&(iscell(x)), ...
'Specify species modification using cell matrix of size [n x 2]');
addParameter(p,'ModSpecies',[],valid_species);
% Optional: simulation options
default_opt = struct;
addParameter(p,'Options',default_opt,@isstruct);
% Optional: parallel vs. series (default is series)
valid_flag = @(x) assert(numel(x)==1, 'Invalid flag passed.');
addParameter(p,'Parallel',1,valid_flag);
% Optional: parallel vs. series (default is series)
addParameter(p,'Plots',1,valid_flag);
% Optional: parallel vs. series (default is series)
expected_stim = {'stim','TAK1'};
addParameter(p,'Stimulus','stim', @(x) any(validatestring(x,expected_stim)));
% Optional: pre-generated simulation data
addParameter(p,'SimData',[], @isnumeric);

% Parse inputs, unpack to variables.
parse(p,doses,varargin{:})
show_graphs = p.Results.Plots;
p_mod = p.Results.ModParams;
species_mod = p.Results.ModSpecies;
stim_name = p.Results.Stimulus;
if strcmp(stim_name,'TAK1')
    p_mod = [p_mod; 66 1 0]; % Clamp TAK1
end

%% SIMULATION: Initialize, run 1st dose (grab steady state), then do remaining simulations (in parallel)

% Specify tracked species and length of simulation (72 hrs)
names = {'IKK','NFkBn'};

if isempty(p.Results.SimData)
    options = p.Results.Options;
    if ~isfield(options,'SIM_TIME')
        options.SIM_TIME = 72*60; % Default simulation time: 72 hrs.
    end
    sim_output = zeros(options.SIM_TIME+1,2, length(doses));

    if show_graphs; disp(['sim #1 - stimulus lvl = ',num2str(doses(1)),'...']); end
    [~,x,simdata] = nfkbSimulate({stim_name,doses(1)},names, p_mod, species_mod,options);
    options.STEADY_STATE = simdata.STEADY_STATE;
    sim_output(:,:,1) = x;

    if p.Results.Parallel
        parfor i = 2:length(doses)
            if show_graphs; disp(['sim #',num2str(i),' - stimulus lvl = ',num2str(doses(i)),'...']); end
            [~,x] = nfkbSimulate({stim_name,doses(i)},names, p_mod, species_mod,options);
            sim_output(:,:,i) = x;
        end
    else
        for i = 2:length(doses)
            if show_graphs; disp(['sim #',num2str(i),' - stimulus lvl = ',num2str(doses(i)),'...']); end
            [~,x] = nfkbSimulate({stim_name,doses(i)},names, p_mod, species_mod,options);
            sim_output(:,:,i) = x;
        end
    end
else
    sim_output = p.Results.SimData;
end

%% FEATURE FINDING
% Find late maximum/minimum (i.e. steady state and oscillation frequency (interpolate to smooth curve)
nfkb_curves = squeeze(sim_output(:,strcmp(names,'NFkBn'),:));
ikk_curves = squeeze(sim_output(:,strcmp(names,'IKK'),:));

maxes = zeros(1,length(doses));
mins = zeros(1,length(doses));
periods = zeros(1,length(doses));
ikk_medians = zeros(1,length(doses));
for i =1:size(nfkb_curves,2)
    maxes(i) = max(nfkb_curves(end-240:end,i));
    mins(i) = min(nfkb_curves(end-240:end,i));
    [~,locs] = findpeaks(nfkb_curves(:,i));
    periods(i) = median(diff(locs));
    ikk_medians(i) = median(ikk_curves(:,i));
end

% Assign feature outputs into single structure
features.steady_state = [mins; maxes];
features.periods = periods;


% Interpolate to smooth curves (for region boundaries, and [optional] graphing)
if length(periods)<2000
    min_d = log10(min(doses));
    max_d = log10(max(doses));
    dose_resamp = 10.^linspace(min_d,max_d,2000);
    maxes = interp1(doses,maxes,dose_resamp);
    mins = interp1(doses,mins,dose_resamp);
    ikk_medians = interp1(doses,ikk_medians,dose_resamp);
    periods = interp1(doses,periods,dose_resamp);
end

% Define region boundaries: pre-oscillatory, oscillations with recovery, damped, and sustained
lo_base = 2; % "Low" baseline: within 2 fold of pre-stimulus value
tol = 5e-3;
bounds(1) = find((maxes-mins)>tol,1,'first');
bounds(2) = find((mins<(lo_base*nfkb_curves(1,1)))& (maxes-mins)>tol,1,'last');
tmp = find((mins>(lo_base*nfkb_curves(1,1)))& (maxes-mins)<tol);
tmp(tmp<bounds(2)) = []; bounds(3) = tmp(1);

%% PLOTS

if show_graphs
    % - - - - GRAPHING PARAMETERS - - - - - - - - - - - - - - - -
    nfkb_lim = [0 0.21];
    dose_lim = [min(doses) max(doses)];
    ikk_lim = [ceil(log10(min(ikk_medians))), (log10(max(ikk_medians)))];
    period_lim = [0 280];
    colors = setcolors;
    grays = colors.bg_gray;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    % #1) Plot 1st 6 hrs of overlaid IKK and NFkB curves
    figs.curves = figure('PaperPositionMode','auto','Position',positionfig(450,112),'Name','IKK and NFkB dynamics');
    ha = tight_subplot(1,2,[0.05 0.05]);
    colororder = cbrewer('div','Spectral',100);
    colororder = colororder([1 5 10 21 71 83 92 100], :);
    set(0,'DefaultAxesColorOrder',colororder)
    dose_show = round(linspace(1,length(doses),8));
    t_show = 0:1/60:6;
    plot(ha(1),t_show,squeeze(sim_output(1:length(t_show),strcmp(names,'IKK'),dose_show)),'LineWidth',2)
    set(ha(1),'xlim',[-0.25 max(t_show)+0.25],'ylim',[-0.005 0.075],'Ytick',0:0.025:0.075,'XTick',0:180:360,...
        'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5)
    grid(ha(1),'on')
    plot(ha(2),t_show,squeeze(sim_output(1:length(t_show),strcmp(names,'NFkBn'),dose_show)),'LineWidth',2)
    set(gca,'xlim',[-0.25 max(t_show)+0.25],'XTickLabel',{},'YTickLabel',{},'YTick',0:0.1:0.3,'XTick',0:180:360,...
        'YLim',nfkb_lim+range(nfkb_lim)*[-0.02 0.3],'LineWidth',1.5)
    grid(ha(2),'on')



    % #2: Feature plots (1) - steady state function of input IKK level
    figs.bifurcation = figure('PaperPositionMode','auto','Position', positionfig(185,155),...
        'Name','Steady states vs IKK level');
    hold(gca,'on')
    area(gca,[ikk_lim(1),ikk_lim(1),log10(ikk_medians([bounds(1) bounds(1)]))],[nfkb_lim, fliplr(nfkb_lim)],...
        'Facecolor',grays(1,:),'EdgeColor','none')
    area(gca,log10(ikk_medians([bounds(1) bounds(1) bounds(2) bounds(2)])), [nfkb_lim, fliplr(nfkb_lim)],...
        'Facecolor',grays(2,:),'EdgeColor','none')   
    area(gca,log10(ikk_medians([bounds(2) bounds(2) bounds(3) bounds(3)])), [nfkb_lim, fliplr(nfkb_lim)],...
        'Facecolor',grays(3,:),'EdgeColor','none')
    plot(gca,log10(ikk_medians),(mins+maxes)/2,'--','Color','k','LineWidth',1.2)
    plot(gca,log10(ikk_medians),[mins; maxes]','Color',colors.red,'LineWidth',1.5)
    hold(gca,'off')
    grid(gca,'on')
    set(gca,'XLim',ikk_lim,'YLim',nfkb_lim,'XTick',-8:2,'YTick',0:0.1:0.3,'XTick',-6:1:0,...
        'XTickLabel',{},'YTickLabel',{},'LineWidth',1.2,'Box','on','Layer','top')

    % #3) Feature plots (2) - median period vs IKK level
    figs.periods = figure('PaperPositionMode','auto','Position', positionfig(185,155),...
        'Name','Period vs IKK level');
    hold(gca,'on')
    area(gca,[ikk_lim(1),ikk_lim(1),log10(ikk_medians([bounds(1) bounds(1)]))],[period_lim, fliplr(period_lim)],...
        'Facecolor',grays(1,:),'EdgeColor','none')
    area(gca,log10(ikk_medians([bounds(1) bounds(1) bounds(2) bounds(2)])), [period_lim, fliplr(period_lim)],...
        'Facecolor',grays(2,:),'EdgeColor','none')   
    area(gca,log10(ikk_medians([bounds(2) bounds(2) bounds(3) bounds(3)])), [period_lim, fliplr(period_lim)],...
        'Facecolor',grays(3,:),'EdgeColor','none')
    plot(gca,log10(ikk_medians),periods,'LineWidth',1.5,'Color',[86 73 158]/255)
    hold(gca,'off')
    grid(gca,'on')
    set(gca,'XLim',ikk_lim,'YLim',period_lim, 'XTickLabel',{},'YTickLabel',{},'YTick',0:60:300,'XTick',-6:1:0,...
        'LineWidth',1.2,'Box','on','Layer','top')

    % #4: Feature plots (3) - steady state vs input level (not IKK)
    figs.stim_vs_ss = figure('PaperPositionMode','auto','Position', positionfig(185,155),'Name',...
        'Steady state vs stimulus level');
    hold(gca,'on')
    hold(gca,'on')
    area(gca,log10([dose_lim(1),dose_lim(1),dose_resamp([bounds(1),bounds(1)])]),[nfkb_lim,fliplr(nfkb_lim)],...
        'FaceColor',grays(1,:),'EdgeColor','none')
    area(gca,log10(dose_resamp([bounds(1) bounds(1) bounds(2) bounds(2)])), [nfkb_lim, fliplr(nfkb_lim)],...
        'FaceColor',grays(2,:),'EdgeColor','none')   
    area(gca,log10(dose_resamp([bounds(2) bounds(2) bounds(3) bounds(3)])), [nfkb_lim, fliplr(nfkb_lim)],...
        'FaceColor',grays(3,:),'EdgeColor','none')
    plot(gca,log10(dose_resamp),(mins+maxes)/2,'--','Color','k','LineWidth',1.5)
    plot(gca,log10(dose_resamp),maxes,'Color',colors.red,'LineWidth',2)
    plot(gca,log10(dose_resamp),mins,'Color',colors.red,'LineWidth',2)
    hold(gca,'off')
    grid(gca,'on')
    set(gca,'XLim',log10(dose_lim),'YLim',nfkb_lim,'XTick',-10:10,'YTick',0:0.1:0.3,...
        'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5,'Box','on','Layer','top')
    
    % #4: Graph of stimulus level vs ikk level
    figs.stim_vs_ikk = figure('PaperPositionMode','auto','Position', positionfig(300,240),'Name',...
        'IKK level vs stimulus level');
    plot(log10(dose_resamp),ikk_medians,'LineWidth',2)
    set(gca,'XLim',log10(dose_lim),'XTick',-10:10,'XTickLabel',{},'YTickLabel',{},'YTick',0:0.02:0.06,'YLim',[0 0.066],...
    'LineWidth',1.5,'Box','on','XGrid','on','YGrid','on','Layer','top')
else
    figs = [];
end

% Switch bounds to an IKK level, not an index.
regimes.ikk = ikk_medians(bounds);
regimes.idx = (bounds/2000)*length(doses);