%% organize sourcedir
% define and create dirs
basedir='C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao';
bidsdir=fullfile(basedir,'BIDS');
cd(bidsdir);
mkdir('sourcedata');
sourcedir=fullfile(bidsdir,'sourcedata');

% create list of subjects in basedir
cd(basedir);
dirs=string(ls);
idx=startsWith(dirs,["HC","PT"]);
subjs=dirs(idx);
subjs=categorical(subjs);
subjs=categories(subjs);

% create a dir per subject in sourcedir
cd(sourcedir);
for i=1:size(subjs,1)
    spm_mkdir(sourcedir,subjs{i});
end

% move raw functional PARRECS from datadir to correct subjectdir under sourcedir
for j=1:size(subjs,1)
    cd(basedir);
    cd(subjs{j});
    cd original_data;
    dirlist1=ls('*.PAR');
    dirlist2=ls('*.REC');
    dirlist=[dirlist1;dirlist2];
    for k=1:size(dirlist,1)
        movefile(dirlist(k,:),fullfile(sourcedir,subjs{j}));
    end
end

% move raw structural PARRECS from datadir to correct subjectdir under
% sourcedir for HC up to HC_021 (rest copied manually - different data
% structure)
for l=1:15
    cd(sourcedir);
    cd(subjs{l});
    mkdir anat;
    cd(basedir);
    cd(subjs{l});
    cd anatomical;
    dirlist3=ls('*.PAR');
    dirlist4=ls('*.REC');
    dirlista=[dirlist3;dirlist4];
    for m=1:size(dirlista,1)
        movefile(dirlista(m,:),fullfile(sourcedir,subjs{l},'anat'));
    end
end

% move raw structural PARRECS from datadir to correct subjectdir under
% sourcedir for PT
idxp=startsWith(subjs,"PT");
pats=subjs(idxp);
for n=1:size(pats,1)
    cd(sourcedir);
    cd(pats{n});
    mkdir anat;
    cd(basedir);
    cd(pats{n});
    cd anatomical;
    dirlist5=ls('*.PAR');
    dirlist6=ls('*.REC');
    dirlistp=[dirlist5;dirlist6];
    for o=1:size(dirlistp,1)
        movefile(dirlistp(o,:),fullfile(sourcedir,pats{n},'anat'));
    end
end

% rename sourcedata subfolders to comply with BIDS
cd (sourcedir);
hc=dir('HC*');
ptfm=dir('PT_FM_*');
ptibs=dir('PT_IBS_*');
ptfmibs=dir('PT_FMIBS_*');
for i = 1:size(hc,1)
    foldername = hc(i).name;
    subnumber = foldername(1,end-1:end);
    newfoldername = strcat('sub-hc',subnumber);
    movefile(foldername,newfoldername); 
end
for i = 1:size(ptfm,1)
    foldername = ptfm(i).name;
    subnumber = foldername(1,end-1:end);
    newfoldername = strcat('sub-ptfm',subnumber);
    movefile(foldername,newfoldername); 
end
for i = 1:size(ptibs,1)
    foldername = ptibs(i).name;
    subnumber = foldername(1,end-1:end);
    newfoldername = strcat('sub-ptibs',subnumber);
    movefile(foldername,newfoldername); 
end
for i = 1:size(ptfmibs,1)
    foldername = ptfmibs(i).name;
    subnumber = foldername(1,end-1:end);
    newfoldername = strcat('sub-ptfmibs',subnumber);
    movefile(foldername,newfoldername); 
end

%% organize rawdir
% make raw directory
cd (bidsdir);
mkdir rawdata;
rawdir=char(fullfile(bidsdir,'rawdata'));

% create list of subjects in sourcedir
cd(sourcedir);
dirs=string(ls);
idx=startsWith(dirs,["sub"]);
source_subs=dirs(idx);
source_subs=categorical(source_subs);
source_subs=categories(source_subs);

% create a dir per subject in rawdir, and create subfolder structure
cd(rawdir);
for i=1:size(source_subs,1)
    spm_mkdir(rawdir,source_subs{i});
    spm_mkdir(rawdir,source_subs{i},'func');
    spm_mkdir(rawdir,source_subs{i},'anat');
end

%% write the dataset description .json file
cd(rawdir);
authors={'Maaike Van Den Houte' 'Katleen Bogaerts' 'Danielle Jongen' 'Huynh Giao Ly' 'Eline Coppens' 'Peter Van Wambeke' 'Koen Schreurs' 'Jan Tack' 'Omer Van den Bergh' 'Lukas Van Oudenhove'};
authors=char(authors);
descr = struct('Name','Brain Mechanisms of Affective Picture Induced Somatic Symptoms in Health and Functional Somatic Syndromes','BIDSVersion','1.4.0','Authors',authors);          
spm_jsonwrite('dataset_description.json',descr, struct('indent','  '));

%% convert functional PARREC to NIFTII and create .json sidecar files
for sub=1:size(source_subs,1)
    % define subjectdirs
    subsourcefuncdir=char(fullfile(sourcedir,source_subs(sub)));
    subrawfuncdir=char(fullfile(rawdir,source_subs(sub),'func'));
    cd(subsourcefuncdir)
    % convert PARREC to .nii and save .nii in rawdir/subject/func folder
    % type help dicm2nii in Matlab command window for more info
    dicm2nii(subsourcefuncdir,subrawfuncdir,1)
end

% note: dicm2nii also converts files in subfolders, hence anatomical
% PARRECS have also been converted already, but live in the
% rawdir/subject/func folder rather than the anat folder, hence we move
% them
cd(rawdir);
for sub=1:size(source_subs,1)
    subrawfuncdir=char(fullfile(rawdir,source_subs(sub),'func'));
    subrawanatdir=char(fullfile(rawdir,source_subs(sub),'anat'));
    cd(subrawfuncdir);
    struct_2_copy = dir('*3DTFE*');
    for struct=1:size(struct_2_copy,1)
        movefile(struct_2_copy(struct).name,subrawanatdir);
    end
end
    
%% rename converted .nii.gz and .json files, and edit .json files where needed
for sub=1:size(source_subs,1)
    subsourcedir=char(fullfile(sourcedir,source_subs(sub)));
    cd(subsourcedir);
    conditions=readtable('order_runs.txt');
    delete('*.mat');
    subrawdir=char(fullfile(rawdir,source_subs(sub)));
    cd(subrawdir);
    anatdir=char(fullfile(subrawdir,'anat'));
    funcdir=char(fullfile(subrawdir,'func'));
        
    % rename .nii files
    % rename .json files, add mandatory taskname to .json files, and
    % anonymize .json files (i.e. clear SeriesDescription)
    cd(anatdir);
    dirlist = dir('*.nii.gz');
        for i = 1:size(dirlist,1)
            filename = dirlist(i).name;
            newfilename = strcat(source_subs(sub),'_T1w.nii.gz');
            newfilename = char(newfilename);
            movefile(fullfile(anatdir,filename),fullfile(anatdir,newfilename)); 
        end
    dirlist = dir('*.json');
        for i = 1:size(dirlist,1)
            filename = dirlist(i).name;
            newfilename = strcat(source_subs(sub),'_T1w.json');
            newfilename = char(newfilename);
                json=spm_jsonread(dirlist(i).name);
                json.SeriesDescription=[];
                spm_jsonwrite(dirlist(i).name,json);
            movefile(fullfile(anatdir,filename),fullfile(anatdir,newfilename)); 
        end  
    cd(funcdir);
    dirlist = dir('*.nii.gz');
        for run = 1:size(dirlist,1)
            filename = dirlist(run).name;
            runname = char(conditions.order(run));
            newfilename = strcat(source_subs(sub),'_task-IAPS_',runname,'_bold.nii.gz');
            newfilename = char(newfilename);
            movefile(fullfile(funcdir,filename),fullfile(funcdir,newfilename)); 
        end
    dirlist = dir('*.json');
        for run = 1:size(dirlist,1)
            filename = dirlist(run).name;
            taskname = char(conditions.order(run));
            newfilename = strcat(source_subs(sub),'_task-IAPS_',taskname,'_bold.json');
            newfilename = char(newfilename);
            n_slices=50;
            TR=3;
            last_slice = TR - TR/n_slices;
            slicetime = linspace(0, last_slice, n_slices);
                json=spm_jsonread(dirlist(run).name);
                json.SeriesDescription=[];
                json.TaskName="IAPS";
                json.SliceTiming=slicetime;
                spm_jsonwrite(dirlist(run).name,json);
            movefile(fullfile(funcdir,filename),fullfile(funcdir,newfilename)); 
        end
    % move dcmHeaders.mat file (output of conversion) back to sourcedir
    % (otherwise conflict when trying to validate BIDS)
    dirlist = dir('dcmHeaders.mat');
        for i = 1:size(dirlist,1)
            filename = dirlist(i).name;
            movefile(fullfile(funcdir,filename),fullfile(subsourcedir,filename)); 
        end
end

%% write events.tsv files - code adapted from NGI12_extract_timing_from_logfiles_emotion.m by Patrick Dupont
%
% DEFINE CONSTANTS
% ----------------
min_time_between_pause_runs   = 145; % in seconds
min_time_between_scoring_runs = 100; % in seconds
min_time_between_fotos_runs   = 100; % in seconds

TR = 3; % in seconds
nr_conditions = 5; % pause, scoring, negative, neutral, positive

itemlist = {
    'Bang.bmp'
    'Beschaamd.bmp'
    'Teneergeslagen.bmp'
    'Vluggeirriteerd.bmp'
    'Zenuwachtig.bmp'
    'Angstig.bmp'
    'Bedroefd.bmp'
    'Gespannen.bmp'
    'Schuldig.bmp'
    'Vijandig.bmp'
    'Duizeligheid.bmp'
    'Misselijkheid.bmp'
    'Snellehartslag.bmp'
    'Spierpijn.bmp'
    'doorademen.bmp'
    'Benauwd.bmp'
    'Bonzen.bmp'
    'Buikpijn.bmp'
    'Hoofdpijn.bmp'
    'Vermoeidheid.bmp'
};

% EXTRACTION OF CONDITIONS, ONSETS, AND DURATIONS FROM LOGFILES
% -------------------------------------------------------------
for sub=1:size(source_subs,1)
    
    % DETERMINE NUMBER OF RUNS IN SUBRAWDIR
    cd(fullfile(sourcedir,source_subs{sub}));
    runs=readtable('order_runs.txt');
    cd(fullfile(rawdir,source_subs{sub},'func'));
    nr_runs=dir('*.nii.gz');
    nr_runs=size(nr_runs,1);
    
    if ~size(runs.order,1) == nr_runs
        error(strcat('number of functional files in rawdir does not match number of runs in order_runs.txt in sourcedir for_', source_subs{sub}));
    end
    
    % DEFINE LOGFILES (ABSOLUTE PATH)
    Plogfile_untitled = spm_select('FPListRec', fullfile(basedir,subjs{sub},'behavioural'), '._untitled.txt');
    Plogfile_pauses   = spm_select('FPListRec', fullfile(basedir,subjs{sub},'behavioural'), '._pauses.txt');
    Plogfile_results  = spm_select('FPListRec', fullfile(basedir,subjs{sub},'behavioural'), '._results.txt');
    
    % READ THE FILE _UNTITLED.TXT (photo blocks)
    fid   = fopen(char(Plogfile_untitled),'r'); 
    data1 = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');
    fclose(fid);
    % blocknumber = column 3
    % foto_tijd   = column 12
    % pos = column 8
    % neg = column 9
    % neu = column 10
    
    tmp_pos = data1(1,8);
    tmp_pos = tmp_pos{1};
    index_headerline1 = [];
        for i = 1:size(tmp_pos,1)
            tmp = char(tmp_pos{i});
            if strcmp(tmp,'pos') == 1
               index_headerline1 = [index_headerline1; i];
            end
        end
    pos = str2num(cell2mat(tmp_pos(index_headerline1(end)+1:end,1)));
    clear tmp_pos
    
    tmp_neg = data1(1,9);
    tmp_neg = tmp_neg{1};
    neg = str2num(cell2mat(tmp_neg(index_headerline1(end)+1:end,1)));
    clear tmp_neg
    
    tmp_neu = data1(1,10);
    tmp_neu = tmp_neu{1};
    neu = str2num(cell2mat(tmp_neu(index_headerline1(end)+1:end,1)));

    tmp_block = data1(1,3);
    tmp_block = tmp_block{1};
    tmp_block = tmp_block(index_headerline1(end)+1:end,1);
    block_number1 = zeros(size(tmp_block));
        for i = 1:size(tmp_block,1)
            block_number1(i) = str2double(cell2mat(tmp_block(i)));
        end

    tmp_foto = data1(1,12);
    tmp_foto = tmp_foto{1};
    tmp_foto = tmp_foto(index_headerline1(end)+1:end,1);
    onset_foto = zeros(size(tmp_foto));
        for i = 1:size(tmp_foto,1)
            onset_foto(i) = str2double(cell2mat(tmp_foto(i)))/1000; % in seconds
        end

    clear data1 index_headerline1

    blocks1    = unique(block_number1);
    nr_blocks1 = length(blocks1);
    
    % READ THE FILE _PAUSES.TXT (pause blocks)
    fid   = fopen(char(Plogfile_pauses),'r'); 
    data2 = textscan(fid,'%s%s%s%s%s%s%s%s','delimiter','\t');
    fclose(fid);
    % blocknumber = column 3
    % pauze_time  = column 8

    tmp_pause = data2(1,8);
    tmp_pause = tmp_pause{1};
    index_headerline2 = [];
    for i = 1:size(tmp_pause,1)
        tmp = char(tmp_pause{i});
        if strcmp(tmp,'pauze_time') == 1
           index_headerline2 = [index_headerline2; i];
        end
    end
    tmp_block = data2(1,3);
    tmp_block = tmp_block{1};
    tmp_block = tmp_block(index_headerline2(end)+1:end,1);
    block_number2 = zeros(size(tmp_block));
    for i = 1:size(tmp_block,1)
        block_number2(i) = str2double(cell2mat(tmp_block(i)));
    end

    tmp_pause = tmp_pause(index_headerline2(end)+1:end,1);
    onset_pause = zeros(size(tmp_pause));
    for i = 1:size(tmp_pause,1)
        onset_pause(i) = str2double(cell2mat(tmp_pause(i)))/1000; % in seconds
    end

    clear data2 index_headerline2
    
    % READ THE FILE _RESULTS.TXT (scoring blocks)
    fid   = fopen(char(Plogfile_results),'r'); 
    data3 = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');
    fclose(fid);
    % blocknumber      = column 3
    % global_start_pos = column 19 (onset in ms)
    % trial contents   = column 7 (item for evaluation)
    % vraagA_1 = column 10
    % vraagB_1 = column 13
    % value    = column 17 (score)

    tmp_vraagA_1 = data3(1,10);
    tmp_vraagA_1 = tmp_vraagA_1{1};
    index_headerline3 = [];
    for i = 1:size(tmp_vraagA_1,1)
        tmp = char(tmp_vraagA_1{i});
        if strcmp(tmp,'vraagA_1') == 1
           index_headerline3 = [index_headerline3; i];
        end
    end
    tmp_vraagA_1 = tmp_vraagA_1(index_headerline3(end)+1:end,1);
    vraagA_1 = zeros(size(tmp_vraagA_1));
    for i = 1:size(tmp_vraagA_1,1)
        vraagA_1(i) = str2double(cell2mat(tmp_vraagA_1(i))); % in seconds
    end

    tmp_block = data3(1,3);
    tmp_block = tmp_block{1};
    tmp_block = tmp_block(index_headerline3(end)+1:end,1);
    block_number3 = zeros(size(tmp_block));
    for i = 1:size(tmp_block,1)
        block_number3(i) = str2double(cell2mat(tmp_block(i)));
    end

    tmp_vraagB_1 = data3(1,13);
    tmp_vraagB_1 = tmp_vraagB_1{1};
    tmp_vraagB_1 = tmp_vraagB_1(index_headerline3(end)+1:end,1);
    vraagB_1 = zeros(size(tmp_vraagB_1));
    for i = 1:size(tmp_vraagB_1,1)
        vraagB_1(i) = str2double(cell2mat(tmp_vraagB_1(i)));
    end

    tmp_score = data3(1,17);
    tmp_score = tmp_score{1};
    tmp_score = tmp_score(index_headerline3(end)+1:end,1);
    score = zeros(size(tmp_score));
    for i = 1:size(tmp_score,1)
        score(i) = str2double(cell2mat(tmp_score(i)));
    end

    tmp_onset = data3(1,19);
    tmp_onset = tmp_onset{1};
    tmp_onset = tmp_onset(index_headerline3(end)+1:end,1);
    onset_scoring = zeros(size(tmp_onset));
    for i = 1:size(tmp_onset,1)
        onset_scoring(i) = str2double(cell2mat(tmp_onset(i)))/1000; % in seconds
    end

    items = data3(1,7);
    items = items{1};
    items = items(index_headerline3(end)+1:end,1);
    clear data3 index_headerline3

    blocks3    = unique(block_number3);
    nr_blocks3 = length(blocks3);
    
    % EXTRACT RELEVANT INFO
    nr_items = length(itemlist);

    % determine the run to which each block belongs
    % rules: onset_foto
    %       1) time between onset_foto is much higher than the time between
    %          two blocks of pictures. We consider values > 100s as differences 
    %          between runs.
    %       2) time between onset_pause is much higher than the time between
    %          two blocks of pictures. We consider values > 145s as differences 
    %          between runs.
    %       3) time between onset_scoring is much higher than the time between
    %          two blocks of scoring. We consider values > 100s as differences 
    %          between runs.
    tmp1 = find(diff(onset_foto) > min_time_between_fotos_runs); % defines the indices of the last part of a run
    tmp2 = find(diff(onset_pause) > min_time_between_pause_runs); % defines the indices of the last part of a run
    tmp3 = find(diff(onset_scoring) > min_time_between_scoring_runs); % defines the indices of the last part of a run

    runs1 = zeros(size(block_number1));
    runs1(1:tmp1(1)) = 1;
    for i = 1:length(tmp1)-1
        runs1(tmp1(i)+1:tmp1(i+1)) = i+1;
    end
    runs1(tmp1(end)+1:end) = i+2;

    runs2 = zeros(size(block_number2));
    runs2(1:tmp2(1)) = 1;
    for i = 1:length(tmp2)-1
        runs2(tmp2(i)+1:tmp2(i+1)) = i+1;
    end
    runs2(tmp2(end)+1:end) = i+2;

    runs3 = zeros(size(block_number3));
    runs3(1:tmp3(1)) = 1;
    for i = 1:length(tmp3)-1
        runs3(tmp3(i)+1:tmp3(i+1)) = i+1;
    end
    runs3(tmp3(end)+1:end) = i+2;
    
    % PUT DATA IN CELL ARRAY FOR WRITING SCORES TO EXCEL
    % NOTE: not all of this is needed for our purpose here, but I kept it
    % in its entirety (except for the actual writing to excel) to avoid
    % breaking the following block of Patrick's code which depends on it -
    % lazy programmer habits, apologies for the "spaghetti code" (term
    % courtesy of Bogdan Petre @CANlab)
    
    output_xls = cell(1+nr_blocks3,25);
    output_xls{1,1} = 'run';
    output_xls{1,2} = 'block number';
    output_xls{1,3} = 'block valence';
    output_xls{1,4} = 'subset_affect';
    output_xls{1,5} = 'subset_somatic';
    for i = 1:nr_items
        clear tmp name
        tmp = char(itemlist(i));
        [~,name,~] = fileparts(tmp);
        output_xls{1,5+i} = name;
    end

    data_scoring = zeros(nr_blocks3,3);
    valence_scoring_pos = zeros(nr_blocks3,1);
    valence_scoring_neg = zeros(nr_blocks3,1);
    valence_scoring_neu = zeros(nr_blocks3,1);
    for i = 1:nr_blocks3
        clear current_block tmp1_indices tmp1_vraagA_1 tmp1_score tmp1_vraagB_1
        clear tmp1_items
        current_block = blocks3(i);
        output_xls{1+i,2} = current_block;
        tmp1_indices = find(block_number3 == current_block);
        tmp1_score = tmp_score(tmp1_indices);
        tmp1_vraagA_1 = vraagA_1(tmp1_indices)+1;
        tmp1_vraagB_1 = vraagB_1(tmp1_indices)+1;
        tmp1_items = items(tmp1_indices);
        if max(unique(tmp1_vraagA_1)) == 2
           output_xls{1+i,4} = 2;
        else
           output_xls{1+i,4} = 1;        
        end
        if max(unique(tmp1_vraagB_1)) == 2
           output_xls{1+i,5} = 2;
        else
           output_xls{1+i,5} = 1;        
        end
        for j = 1:nr_items
            clear current_item tmp2_list tmp2_indices
            current_item = char(itemlist(j));
            tmp2_list = strcmp(tmp1_items,current_item);
            tmp2_indices = find(tmp2_list == 1);
            if length(tmp2_indices) == 1
               output_xls{1+i,5+j} = str2double(cell2mat(tmp1_score(tmp2_indices)));
            end
        end
        data_scoring(i,1) = runs3(tmp1_indices(1));
        data_scoring(i,2) = block_number3(tmp1_indices(1));
        data_scoring(i,3) = onset_scoring(tmp1_indices(1));
        % determine the block valence (neu, pos or neg)
        % find it based on the previous block in blocks1
        indexlist1 = find(block_number1 == current_block-1);
        if pos(indexlist1(1)) == 1
           output_xls{1+i,3} = 'pos';
           valence_scoring_pos(i) = 1;
        elseif neg(indexlist1(1)) == 1
           output_xls{1+i,3} = 'neg';
           valence_scoring_neg(i) = 1;
        elseif neu(indexlist1(1)) == 1
           output_xls{1+i,3} = 'neu';
           valence_scoring_neu(i) = 1;
        end
        output_xls{1+i,1} = unique(runs3(tmp1_indices));
    end
    
    %
    data_foto = zeros(nr_blocks1,3);
    valence_foto_pos   = zeros(nr_blocks1,1);
    valence_foto_neg   = zeros(nr_blocks1,1);
    valence_foto_neu   = zeros(nr_blocks1,1);
    for i = 1:nr_blocks1
        current_block = blocks1(i);
        tmp1_indices = find(block_number1 == current_block);
        data_foto(i,1) = runs1(tmp1_indices(1));
        data_foto(i,2) = block_number1(tmp1_indices(1));
        data_foto(i,3) = onset_foto(tmp1_indices(1));
        if pos(tmp1_indices(1)) == 1
           valence_foto_pos(i)   = 1;
        elseif neg(tmp1_indices(1)) == 1
           valence_foto_neg(i)   = 1;
        elseif neu(tmp1_indices(1)) == 1
           valence_foto_neu(i)   = 1;
        end
    end
    clear current_block current_item
    
    %
    for sess = 1:nr_runs
        % fill in condition onset and duration per run

        % calculate duration of each block within a run
        % assuming that there is always the order: pause foto scoring
        index_foto     = find(data_foto(:,1) == sess);
        index_pause    = find(runs2 == sess);
        index_scoring  = find(data_scoring(:,1) == sess);
        block_foto     = data_foto(index_foto,2);
        block_pause    = block_number2(index_pause);
        block_scoring  = data_scoring(index_scoring,2);
        all_blocks     = [block_pause; block_foto; block_scoring];
        reference_time = onset_pause(index_pause(1));
        onsets         = [onset_pause(index_pause); data_foto(index_foto,3); data_scoring(index_scoring,3)]-reference_time; 
        [order_blocks,order_index] =  sort(all_blocks);
        order_onsets    = onsets(order_index);
        duration_blocks = diff(order_onsets);
        % we don't know exactly the duration of the last scoring period, so I 
        % take the average of all scoring blocks in that run.
        duration_blocks(end+1) = mean(duration_blocks(3:3:end));
        
        % condition pauze
        [~,index_pause,~]  = intersect(order_blocks,block_pause);

        % condition pict_pos
        [~,index_pict_pos,~]  = intersect(order_blocks,blocks1(((data_foto(:,1) == sess).*valence_foto_pos)>0));    

        % condition pict_neu
        [~,index_pict_neu,~]  = intersect(order_blocks,blocks1(((data_foto(:,1) == sess).*valence_foto_neu)>0));    

        % condition pict_neg
        [~,index_pict_neg,~]  = intersect(order_blocks,blocks1(((data_foto(:,1) == sess).*valence_foto_neg)>0));    

        % condition scoring
        [~,index_scoring,~]  = intersect(order_blocks,blocks3(((data_scoring(:,1) == sess))));
        
        trial_type=cell(size(order_onsets,1),1);
        trial_type(index_pause)={'pause'};
        trial_type(index_scoring)={'scoring'};
        trial_type(index_pict_neg)={'negative'};
        trial_type(index_pict_pos)={'positive'};
        trial_type(index_pict_neu)={'neutral'};
        events=table(order_onsets,duration_blocks,trial_type,'Variablenames',{'onset','duration','trial_type'});
        filename=strcat(source_subs{sub},'_task-IAPS_',runs.order{sess},'_events.tsv');
        writetable(events,filename,'Filetype','text','Delimiter','\t');

        % check if all conditions were present in the run, throw a warning
        % if not
        conditionlist = [];
        conditions=unique(events.trial_type);
        for j = 1:nr_conditions
            if ~isempty(conditions(j))
               conditionlist = [conditionlist j];
            else
               warning(['missing condition ' conditions(j) ' in run ' num2str(sess) ])
            end
        end
    end
end