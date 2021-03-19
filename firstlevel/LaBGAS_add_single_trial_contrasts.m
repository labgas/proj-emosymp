%% call function to create DSGN structure array, and define subjectdirs

DSGN = LaBGAS_get_single_trial_dsgn_obj();
sid = dir([DSGN.modeldir,'\sub-*']);
sid = {sid(:).name};
for i = 1:size(sid,2)
    subjmodeldirs{i} = [DSGN.modeldir '\' sid{i}];
end


%% loop over subjects

for j = 2:size(subjmodeldirs,2)
    
    load ([subjmodeldirs{j} '\SPM.mat']);

    betas=SPM.xX.name;
    
    for q = 1:size(DSGN.singletrials{1},2)
        if DSGN.singletrials{1}{q}
            betas_cond{q}=betas(contains(betas,DSGN.conditions{1}{q}));
        end
    end
    
    betas_ofint = {}; % initiate
    r=1;
    
    while r <= size(betas_cond,2)
        betas_ofint=[betas_ofint,betas_cond{r}];
        r=r+1;
    end

    weights=cell(1,size(betas_ofint,2)); % preallocate
    connames=cell(1,size(betas_ofint,2));

    for k = 1:size(betas_ofint,2)
        weights{k}=strcmp(betas_ofint{k},betas);
        weights{k}=double(weights{k});
        connames{k}=betas_ofint{k};
        betas_ofint{k}={{betas_ofint{k}}};
    end

    [matlabbatch, connames, contrast_vectors] = canlab_spm_contrast_job_single_trials_lukasvo(subjmodeldirs{j},betas_ofint,'weights',weights,'names',connames,'nodelete');
    
end


%% subfunction: canlab_spm_contrast_job_single_trials_lukasvo
% function adapted from canlab_spm_contrast_job_luka

function [matlabbatch, connames, contrast_vectors] = canlab_spm_contrast_job_single_trials_lukasvo(modeldir, input_contrasts, varargin)

% set default options
deleteold = 1;
savejob = 1;
runjob = 1;
addcons = 0;
default_suffix = '';

% setup empty cell arrays
custom_connames = cell(numel(input_contrasts),1);
weights = cell(1,numel(input_contrasts));

[matlabbatch, connames, contrast_vectors] = deal({});

% parse varargin
i=1;
while i<=length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'nodelete', deleteold = 0;
            case 'nosave', savejob = 0;
            case 'norun', runjob = 0;
            case 'addcons', addcons = 1; deleteold = 0;
            case 'suffix', i=i+1; default_suffix = varargin{i};
            case 'weights'
                i=i+1;
                for j=1:numel(varargin{i})
                    if ~isempty(varargin{i}{j})
                        weights{j} = varargin{i}{j};
                    end
                end
            case 'names'
                i=i+1;
                for j=1:numel(varargin{i})
                    custom_connames{j} = varargin{i}{j};
                end
            otherwise, error(['Unknown input string option: ' varargin{i}]);
        end
    else
        disp(varargin{i})
        error('Above input argument is unrecognized')
    end
    i=i+1;
end

% setup dir and spm.mat file
spmmatname = fullfile(modeldir,'SPM.mat');

if ~exist(modeldir, 'dir') || ~exist(spmmatname, 'file')
    fprintf('Skipping: No directory or no SPM.mat: %s\n', spmmatname);
    return
end

fprintf('\nContrast Specification for \n%s\n---------------------------------------\n', spmmatname)

% get beta names from SPM.mat file for sanity check
load(spmmatname)
names = SPM.xX.name;
for i = 1:size(input_contrasts,2)
    if ~ismember(input_contrasts{1,i}{1,1},names)
        fprintf('Skipping: Beta name not found in SPM.xX.name: %s\n', input_contrasts{1,i}{1,1}{1,1});
        return
    else
        fprintf('Beta name %s in input_contrasts found in SPM.xX.name - proceeding \n', input_contrasts{1,i}{1,1}{1,1});
    end
end

if addcons
    existing_connames = {};
    for i = 1:numel(SPM.xCon)
        existing_connames{i} = SPM.xCon(i).name; %#ok
    end
end

% setup matlabbatch and attach the SPM mat file
matlabbatch = {};

matlabbatch{1}.spm.stats.con.spmmat(1) = {spmmatname};

% delete previous contrasts if specified
if deleteold    
    matlabbatch{1}.spm.stats.con.delete = 1; % Yes = 1, No = 0;    
else    
    matlabbatch{1}.spm.stats.con.delete = 0; % Yes = 1, No = 0;   
end

% loop over input contrasts
c=0;
for i = 1:size(input_contrasts,2)

    mycon = weights{i};
    conname = custom_connames{i};
    
    if isempty(mycon) % happens if there are missing conditions
        warning('Missing conditions for contrast %s, skipping...', custom_connames{i});
        continue
    end
    
    if addcons
        if numel(existing_connames)>=i & ~strcmp(conname,existing_connames{i})
            warning('input contrast %d (%s) does not match existing contrast %d (%s).',i,conname,i,existing_connames{i})
        end
        
        if strcmp(conname,existing_connames)
            fprintf('Contrast %s exists: skipping.\n',conname)
            continue
        end
    end
    c=c+1;

    if deleteold
        contrastnum = c;
    else
        contrastnum = c + length(SPM.xCon);
    end
    fprintf('Contrast %04d: %s %3.0f Pos weights, %3.0f Neg weights\n', contrastnum, conname, sum(mycon>0), sum(mycon<0));
    
    % Check
    if (any(mycon < 0) && any(mycon > 0)) && (sum(mycon) > 0.00000001)
        disp('WARNING!!!! CONTRAST WEIGHTS DO NOT SUM TO ZERO.  CHECK SPECIFICATION.');
    end
        
    % Save for output
    contrast_vectors{c} = mycon;    
    
    % Add to matlabbatch    
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = conname;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = mycon;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
end % loop over subjects


if c == 0
    fprintf('\nNo contrasts to run: exiting.\n')
    return
end

% save
if savejob    
    savename = fullfile(modeldir,'spm_specify_contrasts_job_single_trials.mat'); % save this as a separate job than the one generated by canlab_spm_contrast_job_luka containing the standard contrasts
    save(savename, 'matlabbatch');
    fprintf('\nSaved batch job (matlabbatch) in spm_specify_contrasts_job_single_trials.mat\n');    
end

% run
if runjob    
    spm('defaults', 'fmri');
    spm_jobman('initcfg'); % initialize
    spm_jobman('run', matlabbatch); % run    
end

end % function