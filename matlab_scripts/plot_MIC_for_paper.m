function [all_diffs,all_baseline] = plot_MIC_for_paper(drug)
    if nargout==1
        all_baseline=0; 
    elseif nargout==2 && ~exist('all_baseline','var')
        all_baseline=0;
    end
    
    close all
    
    % drugs
    drugs = {'Amikacin','Cefepime','Ceftazidime','Ciprofloxacin','Gentamicin','Levofloxacin','Meropenem','Tobramycin'};
    patients = {'13','14','09','42','01','46','10'}; 
    patientlabels = {'B','C','D','F*','G*','H*','I*'};
    
    [S,I,R,minMIC,maxMIC] = get_drug_sensitivity_regimes(drug);
    
    % GET CLINICAL MICs
    clinicalMICs = get_clinical_measurements(); 
    if isfield(clinicalMICs,drug)
        clinicalMICs_drug = clinicalMICs.(drug); 
    else
        clinicalMICs_drug = []; 
    end
    
    OFFSET = 0.1; 
    JITTERSIZE_X = 0.07; 
    JITTERSIZE_Y = 0.09; 
    Y_MIN = minMIC-0.05; 
    Y_MAX = maxMIC+10;
    Y_RANGE = 2.^[-3:6]; 
    
    PLOTTYPE = 'SCATTER';
    fh = figure; hold on; 
    
    if strcmp(PLOTTYPE,'SCATTER')
        samplecolors = [33 146 164; 0 3 128]./255;  
        if strcmp(drug,'KBPipTaz') || strcmp(drug,'KB Cefepime')
            hf2 = fill([0 numel(patients)+1 numel(patients)+1 0], [max(I)*[1 1] Y_MIN*[1 1]], 0.7*[1 1 1]); 
        else
            hf2 = fill([0 numel(patients)+1 numel(patients)+1 0], [I*[1 1] Y_MAX*[1 1]], 0.7*[1 1 1]); 
        end
        set(hf2,'facealpha',.2,'edgecolor','none'); 
        plot([0 numel(patients)+1], minMIC*[1 1],'k--')
        plot([0 numel(patients)+1], maxMIC*[1 1],'k--')
        
    elseif strcmp(PLOTTYPE,'BAR')
        COLORS = [1 1 1; 0 0 0]; 
    end
    
    all_diffs = []; 
    all_baseline = [];
    all_dispersion = [];
    for p = 1:numel(patients)
        [patient_paired_concs,patient_samples] = get_patient_sample_data_for_drug(patients{p},drug);
        if strcmp(PLOTTYPE,'SCATTER')
            numvals_x = size(patient_paired_concs,2); 

            conc_t1 = patient_paired_concs(1,:); 
            plot(get_x_jitter_values(p-OFFSET,numvals_x,JITTERSIZE_X),...
                get_y_jitter_values(conc_t1,JITTERSIZE_Y.*conc_t1),...
                                            'o','markersize',3,'color',samplecolors(1,:),'linewidth',1); 
                                        
            conc_t2 = patient_paired_concs(2,:);
            plot(get_x_jitter_values(p+OFFSET,numvals_x,JITTERSIZE_X),...
                get_y_jitter_values(conc_t2,JITTERSIZE_Y.*conc_t2),...
                                            'o','markersize',3,'color',samplecolors(2,:),'linewidth',1); 
            
            % plot mean within each time point
            barlen = 0.098; 
            plot([p-OFFSET-barlen, p-OFFSET+barlen], nanmean(conc_t1(conc_t1>0))*[1 1], '-',...
                'color',[192 27 27]./255,'linewidth',2);
            plot([p+OFFSET-barlen, p+OFFSET+barlen], nanmean(conc_t2(conc_t2>0))*[1 1], '-',...
                'color',[192 27 27]./255,'linewidth',2);
            
            % mann whitney U test
            [p_val,~]=ranksum(conc_t1(conc_t1>0),conc_t2(conc_t2>0));
            fprintf('Patient %s MW test %0.6f\n',patients{p},p_val);
            
            t1_mean = mean(conc_t1(conc_t1>0)); 
            t2_mean = mean(conc_t2(conc_t2>0)); 
            mean_diff = t2_mean./t1_mean;
            
            all_diffs = [all_diffs mean_diff]; 
            all_baseline = [all_baseline t1_mean]; 
            
            if ~isempty(clinicalMICs_drug)
                try
                    MIC_t1 = clinicalMICs_drug(patient_samples{1}); 
                    MIC_t2 = clinicalMICs_drug(patient_samples{2}); 
                catch 
                    fprintf('No clinical data for %s\n',patients{p});
                end
            end
            
        elseif strcmp(PLOTTYPE,'BAR')
            [n1,x1] = hist(patient_paired_concs(1,:),Y_RANGE);
            [n2,x2] = hist(patient_paired_concs(2,:),Y_RANGE); 
            disp(p*(numel(Y_RANGE)+1) + log2(Y_RANGE)) 
            bh = bar( p*(numel(Y_RANGE)+1) + log2(Y_RANGE), [n1;n2]','grouped'); 
            for b = 1:numel(bh) 
                set(bh(b),'Facecolor',COLORS(b,:))
            end
        end
    end
    
    % prettify plot 
    if strcmp(PLOTTYPE,'SCATTER')
        a1 = gca; 
        set(a1,'xtick',1:numel(patients),'xticklabel',patientlabels,'fontsize',20); 
        
        if strcmp(drug,'KBPipTaz') 
            set(a1,'yscale','linear','ytick',[2 14 21])
            ylim([minMIC maxMIC])
            ylabel('Zone of Inhibition (mm)','fontsize',22);    
        elseif strcmp(drug,'KB Cefepime')
            title('Cefepime','fontsize',22) 
            set(a1,'yscale','linear','ytick',[2 14 18])
            ylim([minMIC maxMIC])
            ylabel('Zone of Inhibition (mm)','fontsize',22);    
        else
            set(a1,'yscale','log','ytick',Y_RANGE); 
            ylim([Y_MIN Y_MAX])
            ylabel(sprintf('%s\nMIC [µg/mL]',drug),'fontsize',22);    
        end
        xlim([0.5 numel(patients)+0.5])
        xlabel('Patients','fontsize',22); 
        
        a2 = axes('YAxisLocation', 'Right','XAxisLocation','bottom');
        set(a2, 'ActivePositionProperty','position')
        set(a2, 'color','none', 'fontsize',20,'xtick',[]); 
        set(a2, 'xtick',[],'xticklabel',[]); 
        box off;
                
        if strcmp(drug,'KBPipTaz') || strcmp(drug,'KB Cefepime')
            set(a2,'ytick',sort([S R],'ascend'),'yticklabel',{'R','S'});
            set(a2,'yscale','linear')
            set(a2,'ylim',[minMIC maxMIC]); 
        else
            set(a2,'ytick',[S R],'yticklabel',{'S','R'}); 
            set(a2,'yscale','log');
            set(a2,'ylim',[Y_MIN Y_MAX]); 
        end
                
        addlistener(a1,'MarkedClean',@(varargin)set(a2,'Position',get(a1,'Position')));
        
    elseif strcmp(PLOTTYPE,'BAR')
        disp('bar') 
    end
    
    
    % --- SUBFUNCTIONS --- %
    
    function xvals = get_x_jitter_values(xval,numpts,jittersize)
        xvals = xval-jittersize + 2*jittersize.*rand(1,numpts);
    end

    function yvals = get_y_jitter_values(yval,jittersize)
        numpts = length(yval); 
        yvals = yval + jittersize.*(-1 + 2*rand(1,numpts));
    end
    
    function [pat_concs,sample_names] = get_patient_sample_data_for_drug(pat,drug)
        datadir = '../MIC_data'; 
        pat_files = dir([datadir filesep pat 'SP*']);
        sample_names = strrep({pat_files.name},'.csv',''); 
        pat_concs = zeros(numel(pat_files),0); 
        for pf = 1:numel(pat_files) 
            alldrugconcs = load_BCH_MIC_sample_data([datadir filesep pat_files(pf).name]); 
            cur_drug = alldrugconcs(drug); 
            pat_concs(pf,1:numel(cur_drug)) = cur_drug'; 
        end
    end

    function drugconcs = load_BCH_MIC_sample_data(fname)
        fid = fopen(fname,'rt'); 
        header = fgetl(fid); 
        drugnames = strsplit(header,',');
        drugconcs = containers.Map(); 
        for d = 1:numel(drugnames)
            drugconcs(drugnames{d}) = []; 
        end
        
        % read sample CSV file 
        while ~feof(fid)
            line = fgetl(fid); 
            comma_placeholdered_line = strrep(line,',,',', ,'); 
            a = strsplit(comma_placeholdered_line,',');
            for col = 2:numel(a)
                % strip extraneous characters
                conc_a = regexprep(a{col},'>|<|\=','');
                drugconcs(drugnames{col}) = [drugconcs(drugnames{col}), str2num(conc_a)];
            end
        end
                
    end

    function clinicalMICs = get_clinical_measurements()
        fid = fopen('../MIC_data/CLINICAL_MICs.csv','rt'); 
        header = fgetl(fid); 
        drugnames = strsplit(header,',');
        clinicalMICs = struct();
        for d = 1:numel(drugnames) 
            clinicalMICs.(drugnames{d}) = containers.Map(); 
        end
        
        % make a struct (drug) of dict (sample name) 
        while ~feof(fid) 
            line = fgetl(fid); 
            comma_placeholdered_line = strrep(line,',,',', ,'); 
            a = strsplit(comma_placeholdered_line,',');
            isolname = a{1}; 
            for col = 2:numel(a)
                % strip extraneous characters
                conc_a = regexprep(a{col},'>|<|\=','');
                if isKey(clinicalMICs.(drugnames{col}),isolname) 
                    clinicalMICs.(drugnames{col})(isolname) = [clinicalMICs.(drugnames{col})(isolname), str2num(conc_a)];
                else
                    clinicalMICs.(drugnames{col})(isolname) = [str2num(conc_a)];
                end 
            end
        end
        
    end

    function [S,I,R,minMIC,maxMIC] = get_drug_sensitivity_regimes(drugname)
        % From CLSI M100 Performance Standards for Antimicrobial Susceptibility Testing 27th edition 
        if strcmp(drugname,'Amikacin')
            S = 16; I = 32; R = 64; 
            minMIC = 2; maxMIC = 64; 
        elseif strcmp(drugname,'Ampicillin')
            minMIC = 1; maxMIC = 32; 
        elseif strcmp(drugname,'Amp/Sul')
            minMIC = nan; maxMIC = 32; 
        elseif strcmp(drugname,'Cefepime')
            S = 8; I = 16; R = 32; 
            minMIC = 1; maxMIC = 64; 
        elseif strcmp(drugname,'Cefoxitin')
            disp('unknown') 
            minMIC = nan; maxMIC = 64; 
        elseif strcmp(drugname,'Ceftazidime')
            S = 8; I = 16; R = 32; 
            minMIC = 1; maxMIC = 64; 
        elseif strcmp(drugname,'Ceftriaxone')
        elseif strcmp(drugname,'Ciprofloxacin')
            S = 1; I = 2; R = 4; 
            minMIC = 0.25; maxMIC = 4; %maxMIC = 64; 
        elseif strcmp(drugname,'Gentamicin')
            S = 4; I = 8; R = 16; 
            minMIC = 1; maxMIC = 64; 
        elseif strcmp(drugname,'Levofloxacin')
            S = 2; I = 4; R = 8; 
            minMIC = 0.125; maxMIC = 64; 
        elseif strcmp(drugname,'Meropenem') 
            S = 2; I = 4; R = 8; 
            minMIC = 0.25; maxMIC = 64; 
        elseif strcmp(drugname,'Tobramycin')
            S = 4; I = 8 ; R = 16; 
            minMIC = 1; maxMIC = 64; 
        elseif strcmp(drugname,'Cefazolin')
            disp('unknown') 
            minMIC = 2; maxMIC = 64; 
        elseif strcmp(drugname,'KB Cefepime')
            S = 18; I = 15:17; R = 14; 
            minMIC = 2; maxMIC = 30; 
        elseif strcmp(drugname,'KBPipTaz')
            S = 21; I = 15:20; R = 14; 
            minMIC = 2; maxMIC = 30; 
        end
    end
    
end