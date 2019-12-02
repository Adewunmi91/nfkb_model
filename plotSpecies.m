function [fig1] =  plotSpecies(simdata, names, doses, t)
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [figs] =  plotSpecies(simdata, t)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PLOTSPECIES takes a simulation run (output of nfkbSimulate) and plots each output species in a subplot
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default t vector: assume 1-minute spacing
if nargin<4
    t = 0:(size(simdata,1)-1);
end

all_color = linspecer(100);
color_order = all_color(round(linspace(1,100, length(doses))),:); %colormap_byr(round(linspace(1,64,length(doses))),:);
n_rows = ceil((length(names)+1)/5);
n_cols = min([5, 1+length(names)]);
t = t/60;
t_pad = range(t)*0.02;

fig1 = figure('Position',positionfig(1000, n_rows*120 + 50),'ResizeFcn',@fix_ticks);
ha = tight_subplot(n_rows,n_cols,[0.08 0.04],[0.08 0.08],[0.04 0.04]);
for i = 1:length(names)
    hold(ha(i),'on')
    set(ha(i),'ColorOrder',color_order)
    vects = squeeze(simdata(:,i,:))';
    plot(ha(i),t,vects,'LineWidth',2)
    min_y = real(max([min(vects(:))-range(vects(:))*0.08, -range(vects(:))*0.03]));
    max_y = real(max(vects(:))+ range(vects(:))*0.1);
    if abs(min_y-max_y)<eps
        max_y = max_y + max([0.02*min_y,eps]);
    end
    hold(ha(i),'off')
    set(ha(i),'XLim',[-t_pad max(t)+t_pad],'FontSize',9,'Box','on','Layer','top','YLim', sort([min_y max_y]))
    labels = get(ha(i),'YTick');
    set(ha(i),'YTickLabel',labels);

    title(ha(i),names{i},'Interpreter','none');
    if i>(length(ha)-5)
        set(ha(i),'XTickLabel',get(ha(i),'XTick'))
    end
    grid(ha(i),'on')
end
labels = cellfun(@(x) num2str(x,'%.2g'),num2cell(doses),'UniformOutput',false);
colormap(color_order)
hcb=colorbar('Ticks',(1:length(doses))+0.5,'TickLabels',labels,'FontSize',10);
hcb.Label.String = 'Stimulus level';
hcb.Label.FontSize = 14;
caxis([1 length(doses)+1])

function fix_ticks(varargin)
ha = get(varargin{1},'Children');
if ~isempty(ha)
    % Axes are reversed from 'tight_subplot' order
    for i = 3:6 % Skip colorbar and empty axes
        set(ha(i),'XTickLabel',get(ha(i),'XTick'),'FontSize',9)
    end
end



