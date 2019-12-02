function [cell_err, cell_params] = scan_to_cell(model_depth)

% 1st 64 rows of (ordered) model parameters (from supp_modelfit)
best_params = [...
4         1         14        1         9         12        
4         1         11        4         9         12        
0.25      1         11        1         6         12        
1         4         14        1         9         12        
1         4         11        4         9         12        
0.25      1         11        1         9         12        
4         4         11        4         6         16        
4         4         14        1         6         16        
1         1         11        1         9         12        
1         0.25      14        1         9         12        
1         0.25      11        4         9         12        
1         1         11        1         6         16        
4         4         14        1         9         12        
4         4         11        4         9         12        
4         1         14        1         9         16        
4         1         11        4         9         16        
4         1         11        4         6         16        
4         1         14        1         6         16        
4         0.25      14        1         9         12        
4         0.25      11        4         9         12        
4         4         11        4         6         12        
4         4         14        1         6         12        
4         0.25      14        4         6         12        
4         0.25      11        7         6         12        
4         0.25      17        1         6         12        
0.25      1         11        1         3         12        
4         4         14        1         9         16        
4         4         11        4         9         16        
4         0.25      11        4         9         16        
4         0.25      14        1         9         16        
1         0.25      11        4         6         12        
1         0.25      14        1         6         12        
1         1         11        4         3         16        
1         1         14        1         3         16        
4         1         17        1         3         12        
4         1         11        7         3         12        
4         1         14        4         3         12        
0.25      4         11        1         9         12        
0.25      4         11        1         6         12        
4         1         17        1         3         16        
4         1         11        7         3         16        
4         1         14        4         3         16        
1         4         14        1         6         16        
1         4         11        4         6         16        
1         1         14        1         3         12        
1         1         11        4         3         12        
1         1         14        1         6         12        
1         1         11        4         6         12        
1         1         11        4         9         12        
1         1         14        1         9         12        
0.25      4         11        1         3         12        
1         4         11        4         6         12        
1         4         14        1         6         12        
1         1         11        1         9         16        
4         4         14        1         3         20        
4         4         11        4         3         20        
4         1         14        1         6         12        
4         1         11        4         6         12        
1         4         14        1         9         16        
1         4         11        4         9         16        
1         4         11        1         9         12        
0.25      1         11        1         3         16        
4         1         14        4         6         12        
4         1         17        1         6         12        
];

nfkb = loadnfkb;
load('tnf_peaks.mat')

% Cells to attempt to match
nfkb_num = [2,3,4,5,6]; % Original TNF nfkb sets that correspond to peak_times, etc.
fit_cells = [
    3   86
    4   174
    4   135
    4   233
    ];


% Create arrays based on number of (ranked) models to evaluate
cell_err = zeros(size(fit_cells,1),model_depth);
cell_params = zeros(5,size(fit_cells,1),model_depth);

for i = 1:size(fit_cells,1)
    % Grab trajectory to match (starting at zero point); interpolate it to minute intervals
    trajectory = medfilt1(nfkb(nfkb_num(fit_cells(i,1))).metrics.time_series(fit_cells(i,2),2:end),4);
    t_max = 335;
    t = 0:5:t_max;
    trajectory = trajectory(1:length(t));
    trajectory = interp1(t,trajectory,0:t_max);
    
    % Attempt to optimize "guess" parameters (IKK-IkBa binding, IkBa induction level) to fit trajectories
    options = psoptimset('MaxFunEvals',250,'Display','final');
    params = [
    %  idx1 idx2  x0    ub      lb
        6   2   2       4     1
        6   3   0.24    1       0.01
        21  1   150     1e3     1e-1
        23  1   8       100     1e-1
        25  1   2       20      1e-1
    ];
    param_idx = params(:,1:2);
    x0 = params(:,3)';
    ub = params(:,4)';
    lb = params(:,5)';
    parfor j = 1:model_depth
        extra_params = [
         % IkBa transport: rxns #9 and #11
            9   1   0.09*best_params(j,2)
            11  1   0.09*best_params(j,2)/2*3.5
            % NFkB transport: rxns #10 and #12
            10  1    0.6*best_params(j,1)
            12  1    0.6*best_params(j,1)/50*3.5
            % discrete delay 1 (mRNA)
            6   4   best_params(j,3)
            % discrete delay 2 (protein)
            8   2   best_params(j,4)
            % IkBa halflife
            15  1   log(2)/best_params(j,5)
            16  1   log(2)/best_params(j,5)
            % IkBat halflife
            7   1   log(2)/best_params(j,6)
            ];
        fitfcn = @(x) fit_oscillations(x,param_idx,trajectory, extra_params);
        tic
        [cell_params(:,i,j),cell_err(i,j)] = patternsearch(fitfcn, x0,[],[],[],[],lb,ub,[],options);
        toc
    end
end


