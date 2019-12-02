% Define experimantal data (Werner 2005, Science)
ikk = struct;
ikk.mef = [1 59 100 64 50 36 21 20 18 16 14 12 10 8 7 5 4 2.5 1 1 1]; % MEF data
ikk.fibro = [1 50 100 30 18 6 5.25 4.5 3.73,3,2.67, 2.33, 2, 1.67 1.33 1 1 1 1 1 1]; % 3T3 data
t_exp = 5*((1:length(ikk.mef))-1);
fit.t = [0:5:50, 90];
fit.ikk = ikk.mef(ismember(t_exp,fit.t));

noise_add = @(x) x*(1000^(rand(1)-0.5));
n_runs = 1000;
all_x = zeros(11,n_runs);
all_err = zeros(1,n_runs);

parfor i = 1:n_runs
    t0 = 0:max(t_exp);
    options = psoptimset('MaxFunEvals',500,'Display','final');
    params = [
    %  idx1 idx2  x0    ub      lb
        58  1   0.08  0.14    0.035
        59  1   10      1000    1
        60  1   0.09    1000    1e-3
        62  1   5       1000    1
        63  1   0.02    1000    1e-4
        65  1   500     2000    10
        66  1   1       100     1e-3
        67  1   1.4     100      0.1
        67  3   0.004   5e-2    5e-4
        3   1   2       10      1e-1
        4   1   10      200     1e-1
    ];

    % Add log-normal noise  w/ high c.v.
    for j = 1:size(params,1)
        params(j,3) = noise_add(params(j,3));
    end
    
    % Cap parameters at lower/upper bounds
    cap_row = params(:,3) > params(:,4);
    params(cap_row,3) = params(cap_row,4);
    cap_row = params(:,3) < params(:,5);
    params(cap_row,3) = params(cap_row,5);
    
    disp('Starting x=')
    spell(params(:,1:3));

    param_idx = params(:,1:2);
    x0 = params(:,3)';
    ub = params(:,4)';
    lb = params(:,5)';
    fitfcn = @(x) fit_ikk_tnf(x,param_idx,fit);
    tic
    [p_fit,fval] = patternsearch(fitfcn, x0,[],[],[],[],lb,ub,[],options);
    all_err(i) = fval;
    all_x(:,i) = p_fit';
toc
end

save('/home/brooks/Data/ikk_scan.mat','all_err','all_x')

%%

params = [
%  idx1 idx2  x0    ub      lb
    58  1   0.1  0.14    0.035
    59  1   10      1000    1
    60  1   0.09    1000    1e-3
    62  1   5       1000    1
    63  1   0.02    1000    1e-4
    65  1   500     2000    10
    66  1   1       100     1e-3
    67  1   1.4     100      0.1
    67  3   0.004   5e-2    5e-4
    3   1   2       10      1e-1
    4   1   10      200     1e-1
];

param_idx = params(:,1:2);

% Define experimantal data (Werner 2005, Science)
ikk = struct;
ikk.mef = [1 59 100 64 50 36 21 20 18 16 14 12 10 8 7 5 4 2.5 1 1 1]; % MEF data
ikk.fibro = [1 50 100 30 18 6 5.25 4.5 3.73,3,2.67, 2.33, 2, 1.67 1.33 1 1 1 1 1 1]; % 3T3 data
t_exp = 5*((1:length(ikk.mef))-1);
fit.t = [0:5:50, 90];
fit.ikk = ikk.mef(ismember(t_exp,fit.t));

depth = 10;

[~,order] = sort(all_err,'ascend');
for i = 1:depth
    p_fit = all_x(:,order(i))';
    [err,ikk] = fit_ikk_tnf(p_fit,param_idx,fit,1);
    disp(['Param set #',num2str(i)])
    spell(p_fit')
    disp(['error = ',num2str(err)])
    all_ikk(i) = ikk(46);
    disp('- - - - ')
end