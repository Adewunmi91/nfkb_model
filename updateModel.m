function [] = updateModel(spreadsheet)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = updateModel(spreadsheet)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% UPDATEMODEL rewrites default parameters and ODE equations from the "whole-cell NFkB model" spreadsheet
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%
% Google spreadsheet location

if nargin<1
    spreadsheet = 'https://docs.google.com/spreadsheets/d/1r1YOpwae7e93TAHzQwaYj-hR1VasCU40h2jxBLr4yWw/pubhtml';
    disp('Reading model from default Google spreadsheet:')
    disp(spreadsheet)
end

if exist(spreadsheet,'file')
    [rxn_data, species_data, param_num, param_str, col_start] = parseExcel(spreadsheet);
else
    [rxn_data, species_data, param_num, param_str, col_start] = parseGoogle(spreadsheet);
end


%% Form reaction equations  (save any volume modifiers for "delta" equations later)
rxn_eqns = {};
all_vol ={};
for i =1:size(param_num,1)
    if param_num(i,2)==1
        reactants = strsplit(rxn_data{i,col_start(i)},'+');
        rxn_idx = num2str(param_num(i,1));
        vol_scale = {'',''};
        % Parse reaction modifiers
        modifier_set = param_str(param_num(:,1)==param_num(i,1),1);
        modifier_set(1) = [];
        % MODIFIER #1: delay -> if found, modify ALL reactant species
        if min(cellfun(@isempty, strfind(lower(modifier_set),'delay')))==0
            for k = 1:length(reactants)
                if ~isempty(strtrim(reactants{k}))
                    reactants{k} = [strtrim(reactants{k}),'_tau',num2str(param_num(i,1))];
                end
            end
        end
        % MODIFIER #2: Hill -> if found, modify FIRST reactant species
        hill_idx = find(cellfun(@isempty, strfind(lower(modifier_set),'hill'))==0,1,'first');
        if ~isempty(hill_idx)
            if ~isempty(strtrim(reactants{1}))
                hill_idx = num2str(hill_idx+1);
                km_idx = num2str(find(cellfun(@isempty, strfind(lower(modifier_set),'km'))==0,1,'first')+1);
                x = strtrim(reactants{1});
                reactants{1} = ['(',x,'.^p(',rxn_idx,',',hill_idx,'))/',...
                    '( (',x,'.^p(',rxn_idx,',',hill_idx,')) + '...
                    '(p(',rxn_idx,',',km_idx,').^p(',rxn_idx,',',hill_idx,')) )'];
            end
        end
        % MODIFIER #3: volume scale ->  if found, set  correct row of volume scale
        vol_idx = find(cellfun(@isempty, strfind(lower(modifier_set),'vol.'))==0,1,'first');
        if ~isempty(vol_idx)
            if isempty(strfind(modifier_set{vol_idx},'reactant'))
                vol_scale{2} = ['p(',rxn_idx,',',num2str(vol_idx+1),')'];
            else
                vol_scale{1} = ['p(',rxn_idx,',',num2str(vol_idx+1),')'];
            end
        end
        
        % Assemble full reaction equation 
        rxn_eqn = ['rxn_',rxn_idx,' = p(',rxn_idx,',1)'];
        for j = 1:length(reactants)
            if ~isempty(strtrim(reactants{j}))
                rxn_eqn = [rxn_eqn,' * ',strtrim(reactants{j})];
            end
        end
            rxn_eqns = cat(1,rxn_eqns,[rxn_eqn,';']);
            all_vol = cat(1,all_vol,vol_scale);
    end
end

%% Form delta equations
reactants = {}; products = {};
for i = 1:size(param_num,1)
    if param_num(i,2)==1
        reactants = cat(1,reactants,{strsplit(rxn_data{i,col_start(i)},'+')});
        products = cat(1,products,{strsplit(rxn_data{i,col_start(i)+2},'+')});
    end
end
species_names = species_data(2:end,2);

delta_eqn = cell(size(species_names));
for i = 1:length(delta_eqn)
    % Collect all reactions where each species participates
    r_nums = [];
    p_nums = [];
    r_scale = [];
    p_scale = [];
    for j = 1:length(reactants)
        flag = 0;
        for k=1:length(reactants{j})
            if strcmp(strtrim(reactants{j}{k}),strtrim(species_names{i}))
                r_nums = [r_nums,j];
                if k==1
                    r_scale = [r_scale,1];
                else
                    r_scale = [r_scale,0];
                end
                flag = 1;
            end
        end
        for m=1:length(products{j})
            if strcmp(strtrim(products{j}{m}),strtrim(species_names{i}))
                p_nums = [p_nums,j];
                if m==1
                    p_scale = [p_scale,1];
                else
                    p_scale = [p_scale,0];
                end
                if flag; flag=2; end
            end
        end
        if flag>1
            r_nums = r_nums(1:end-1);
            p_nums = p_nums(1:end-1);
            r_scale = r_scale(1:end-1);
            p_scale = p_scale(1:end-1);
        end       
    end
    % Combine reactions into a single expression
    delta_eqn{i} = ['delta(',num2str(i),') ='];
    if isempty(r_nums) && isempty(p_nums)
        delta_eqn{i} = [delta_eqn{i},' 0'];
    else
        for j = 1:length(r_nums)
            if ~isempty(all_vol{r_nums(j),1}) && r_scale(j)
                addendum = [' - (rxn_',num2str(r_nums(j)),' * ',all_vol{r_nums(j),1},')'];
            else
                addendum = [' - rxn_',num2str(r_nums(j)),''];
            end
            delta_eqn{i} = [delta_eqn{i}, addendum];
        end
        for j = 1:length(p_nums) 
            if ~isempty(all_vol{p_nums(j),2}) && p_scale(j)
                addendum = [' + (rxn_',num2str(p_nums(j)),' * ',all_vol{p_nums(j),2},')'];
            else
                addendum = [' + rxn_',num2str(p_nums(j))];
            end
            delta_eqn{i} = [delta_eqn{i}, addendum];
        end
    end
    delta_eqn{i} = [delta_eqn{i},';'];
end




%% Write files. 
parent_path = mfilename('fullpath');
seps = strfind(parent_path,filesep);
parent_path = parent_path(1:seps(end));

% 1) 'nfkbInitialize' contains physical constants, parameters, and species.
fileID = fopen([parent_path,'nfkbInitialize.m'],'w');
fprintf(fileID,'function [params, species] = nfkbInitialize()\n');
fprintf(fileID,'%% This file is automatically generated by updateModel.m from the NFkB whole-cell spreadsheet\n');
fprintf(fileID,['%% Model URL: ', spreadsheet,'\n\n\n']);
fprintf(fileID,'%% PARAMETERS\n');
formatSpec = 'params(%d,%d) = %g; %% %s\n';
for i = 1:size(param_num,1)
    fprintf(fileID,formatSpec,param_num(i,1),param_num(i,2),param_num(i,3),param_str{i,2});
end
fprintf(fileID,'\n\n%% SPECIES\n');
fprintf(fileID,'species.NAMES = {...\n');
formatSpec = '''%s'' ...%d\n';
for i = 1:length(species_names)
    fprintf(fileID,formatSpec,species_names{i},i);
end
fprintf(fileID,'};');
fclose(fileID);


%% 2) 'nfkbOde' sets up reactions and calculates differentials - pulls in 2 subfiles, ODE_INIT and ODE_DELAY
fileID = fopen([parent_path,'nfkbOde.m'],'w');
fprintf(fileID,'%% This file is automatically generated by updateModel.m from the NFkB whole-cell spreadsheet.\n');
fprintf(fileID,['%% Model URL: ', spreadsheet,'\n']);

% Section 1: Get initialization code from ODE_INIT
fprintf(fileID,'\n%%%% Section 1: Declaration/initialization (code from ode_init.m)\n');
subID = fopen([parent_path,'ode_init.m'],'r');
content = fread(subID, '*char');
fclose(subID);
fwrite(fileID,content);

% Section 2: Get species names, fill in values from 'x' vector
fprintf(fileID,'\n\n\n%%%% Section 2: Unpack species\n');
formatSpec = '%s = x(%d);\n';
for i = 1:length(species_names)
    fprintf(fileID,formatSpec,species_names{i},i);
end

% Section 3: Get delay code from ODE_DELAY
fprintf(fileID,'\n\n%%%% Section 3: Calculate delays (code from ode_delay.m)\n');
subID = fopen([parent_path,'ode_delay.m'],'r');
content = fread(subID, '*char');
fclose(subID);
fwrite(fileID,content);
fprintf(fileID,'\n');

% Section 4: Reaction fluxes
fprintf(fileID,'\n\n%%%% Section 4: Set reaction rates\n');
formatSpec = '%s\n';
for i = 1:size(rxn_eqns,1)
    fprintf(fileID,formatSpec,rxn_eqns{i});
end

% Section 5: Species deltas
fprintf(fileID,'\n\n%%%% Section 5: Set species'' deltas from reactions\n');
formatSpec = '%s\n';
for i = 1:size(delta_eqn,1)
    fprintf(fileID,formatSpec,delta_eqn{i});
end
fclose(fileID);




function [rxn_data, species_data, param_num, param_str, col_start] = parseGoogle(url)
%% GOOGLE SPREADSHEET specified: parse HTML tables
pageString = urlread(url);
[tables] = regexp(pageString, '(<table[^>]*>(?:(?>[^<]+)|<(?!table[^>]*>))*?</table>)','tokens');
titles =  regexp(pageString, '(<li id="sheet-button[^>]*>(?:(?>[^<]+)|<(?!li[^>]*>))*?</li>)','tokens');
% Build cell aray of table data
for i = 1:length(tables)
    table = tables{i};
    rows = regexpi(table, '<tr.*?>(.*?)</tr>','tokens');
    table_data = cell(0);
    % Pull off HTML <th> headers (if present)
    headers = regexpi(rows{1}{1}, '<th.*?>(.*?)</th>','tokens');
    if isempty(headers{1})
        start_mod = 0;
    else
        start_mod = 1;
    end
    % Cycle rows, then columns - pull information into a cell array
    for j = 1:(numel(rows{1})-start_mod)
        cols = regexpi(rows{1}{j+start_mod}, '<td.*?>(.*?)</td>','tokens');
        for k = 1:numel(cols{1})
            tmp = regexprep(cols{1}{k},'<.*?>', '');
            table_data{j,k} = tmp{1};
        end
    end

    % Peel off reactions and species tables - break loop when both are found
    title = titles{i};
    % 'Species' table
    if ~isempty(strfind(title{1},'Species<'))
        species_data = table_data;
    end
    % 'Reactions' table
    if ~isempty(strfind(title{1},'Reactions<'))
        rxn_data = table_data;
    end
    if exist('species_data','var')&&exist('rxn_data','var')
        break
    end
end

% Drop any empty rows/columns
rxn_data(:,sum(cellfun(@isempty,rxn_data),1)==size(rxn_data,1)) = [];
rxn_data(sum(cellfun(@isempty,rxn_data),2)==size(rxn_data,2),:) = [];

species_data(:,sum(cellfun(@isempty,species_data),1)==size(species_data,1)) = [];
species_data(sum(cellfun(@isempty,species_data),2)==size(species_data,2),:) = [];

% Pull hearders off data; identify columns for parameter index and value
rxn_headers = rxn_data(1,:);
rxn_data = rxn_data(2:end,:);
id_col = find(strcmpi('#',rxn_headers),1,'first');
value_col = find(strcmpi('value',rxn_headers),1,'first');

% Collect parameters, values, units, and descriptions. 
% (Merge cell behavior for Google: the first cell is kept, and the others are deleted)
param_num = nan(size(rxn_data,1),3);
param_str = cell(size(rxn_data,1),2);
col_start = 2*ones(size(rxn_data,1),1);

for i = 1:size(param_num,1)
    % Fill in parameter indicies
    try
        param_num(i,1) = eval(rxn_data{i,id_col});
    catch me
        param_num(i,1)  = param_num(i-1,1);
        col_start(i) = 1;
    end
    param_num(i,2) = sum(param_num(:,1)==param_num(i,1));
    % Fill in parameter values
    for j = col_start(i):value_col
        try
            param_num(i,3) = rxn_data{i,j};
            unit_start = j+1;
            break;
        catch me
        end
    end
    % Fill in unit and description
    param_str{i,1} = rxn_data{i,unit_start};
    param_str{i,2} = rxn_data{i,unit_start+1};
end


function [rxn_data, species_data, param_num, param_str, col_start] = parseExcel(spreadsheet)
%% EXCEL SHEET specified: read in reaction data and species data
[~, titles] = xlsfinfo(spreadsheet);
tables = cell(size(titles));
for i = 1:length(titles)
    [~,~,tables{i}] = xlsread(spreadsheet,titles{i});
    % 'Species' table
    if strcmpi(titles{i},'Species')
        species_data = tables{i};
    end
    % 'Reactions' table
    if strcmpi(titles{i},'Reactions')
        rxn_data = tables{i};
    end
    if exist('species_data','var')&&exist('rxn_data','var')
        break
    end
end

% Drop empty rows
numerics = cellfun(@isscalar,cellfun(@isnan,rxn_data,'UniformOutput',false));
numeric_rows = find(min(numerics,[],2) > 0);
rxn_data(numeric_rows(min(cellfun(@isnan,rxn_data(numeric_rows,:)),[],2)>0),:) = [];


numerics = cellfun(@isscalar,cellfun(@isnan,species_data,'UniformOutput',false));
numeric_rows = find(min(numerics,[],2) > 0);
species_data(numeric_rows(min(cellfun(@isnan,species_data(numeric_rows,:)),[],2)>0),:) = [];


% Initialize structures to hold parameter values
rxn_headers = rxn_data(1,:);
rxn_data = rxn_data(2:end,:);
id_col = find(strcmpi('#',rxn_headers),1,'first');
value_col = find(strcmpi('value',rxn_headers),1,'first');

    
 % Collect parameters, values, units, and descriptions. 
% (Merge cell behavior for Excel: the first cell is kept, and the others are converted to NaN)
param_num = nan(size(rxn_data,1),3);
param_str = cell(size(rxn_data,1),2);
col_start = 2*ones(size(rxn_data,1),1);

for i = 1:size(param_num,1)
    % Fill in parameter indicies
    param_num(i,1) = rxn_data{i,id_col};
    if isnan(param_num(i,1))
        param_num(i,1)  = param_num(i-1,1);
    end
    param_num(i,2) = sum(param_num(:,1)==param_num(i,1));
    % Fill in parameter values
    val = rxn_data{i,value_col};
    if isnumeric(val)
        param_num(i,3) = val;
    else
        try
            param_num(i,3) = eval(val);
        catch me
            error(['Bad value passed from spreadsheet. (row ',num2str(i),', row = ''', rxn_data{i,value_col},''')'])
        end
    end
    % Fill in unit and description
    param_str{i,1} = rxn_data{i,value_col+1};
    param_str{i,2} = rxn_data{i,value_col+2};
end


% Convert all rxn_data numbers to empty strings
c = cell(size(rxn_data));
for i = 1:numel(c)
    c{i} = '';
end
rxn_data(cellfun(@isnumeric,rxn_data)) = c(cellfun(@isnumeric,rxn_data));
