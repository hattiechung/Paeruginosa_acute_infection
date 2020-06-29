function plot_treatment_MIC_correlation(drug_class)
    % 'drug_class': can be 'blactam' or 'quinolone'
    %%
    % get patient treatment and MICs for patients with paired samples

    patient_nums = [13 14 9 42 1 46 10]; 
    patient_str = {'13','14','09','42','01','46','10'};  
    patient_letters = {'B','C','D','F*','G*','H*','I*'};
    n_pats = numel(patient_nums); 

    date_span = 0:15; 
    blactam_abx = zeros(numel(patient_nums),numel(date_span));
    cefep_abx = zeros(numel(patient_nums),numel(date_span));
    mero_abx = zeros(numel(patient_nums),numel(date_span));
    cipro_abx = zeros(numel(patient_nums),numel(date_span));
    azithro_abx = zeros(numel(patient_nums),numel(date_span));

    % returned drugs, in order {'cefep','ceftaz','ceftriax','mero','piptaz','cipro','azithro'};
    for p = 1:n_pats
        [patient_med_dates, days_fr_sample1] = plot_MIC_abx_course(patient_nums(p),false);
        if numel(patient_med_dates)>3
            for d = 1:numel(date_span) 
                cur_d = date_span(d); 
                if find(days_fr_sample1==cur_d)
                    blactam_abx(p,d) = sum(patient_med_dates([1:3,5],find(days_fr_sample1==cur_d))); 
                    cefep_abx(p,d) = patient_med_dates(1,find(days_fr_sample1==cur_d));
                    mero_abx(p,d) = patient_med_dates(4,find(days_fr_sample1==cur_d));
                    cipro_abx(p,d) = patient_med_dates(6,find(days_fr_sample1==cur_d));
                    azithro_abx(p,d) = patient_med_dates(7,find(days_fr_sample1==cur_d));
                end
            end
        end
    end

    patient_samples = []; 
    for p = 1:n_pats
        pat = patient_str{p};
        [psamples, sampledays, stypes] = get_patient_sample_dates(pat); 
        patient_samples(p,sampledays) = stypes; 
    end
    
    %% determine which drug to plot
    if strcmp(drug_class,'blactam')
        abx_course = blactam_abx;
        mic_drug = 'Cefepime'; 
    elseif strcmp(drug_class,'quinolone')
        abx_course = cipro_abx;
        mic_drug = 'Ciprofloxacin';     
    end
    
    %% calculate corr between samples
    day1 = find(date_span==0); 
    abx_bt_samples = logical(abx_course(:,day1:day1+size(patient_samples,2)-1));
    treated_days_frac = []; 
    for p = 1:n_pats
        t1 = find(patient_samples(p,:)==1);
        t2 = find(patient_samples(p,:)==2);
        days_treated = nansum(abx_bt_samples(p,t1:t2));
        days_total = numel(t1:t2);
        fprintf('Patient %s %0.3f\n',patient_str{p},days_treated./days_total);
        treated_days_frac = [treated_days_frac days_treated./days_total];
    end

    [percent_diffs, ~] = plot_MIC_for_paper(mic_drug);

    %% plot
    fh_corr = figure; 
    scatter(treated_days_frac, percent_diffs, 30,'markerfacecolor',0.2*[1 1 1],'markeredgecolor','k');
    [R_coeff,p_corr]=corrcoef(treated_days_frac,percent_diffs); 
    box off; 

    % fit line
    P = polyfit(treated_days_frac,percent_diffs,1);
    yfit = P(1)*treated_days_frac+P(2);
    hold on;
    plot(treated_days_frac,yfit,'color',0.2*[1 1 1],'linewidth',2);
    plot([-1 2], [1 1], 'k:', 'linewidth', 1);
    
    for a = 1:numel(percent_diffs)
        text(treated_days_frac(a)+.025,percent_diffs(a)+.1,patient_letters{a},'fontsize',16); 
    end
    
    text(0,2.8,sprintf('R^2 = %0.3f\np = %0.3f',R_coeff(1,2),p_corr(1,2)),'fontsize',18);

    % prettify
    xlabel(sprintf('%% days between samples\ntreated with drug class'),'fontsize',18);
    ylabel(sprintf('Mean fold change\n%s resistance', mic_drug),'fontsize',18);
    xlim([-0.1 1.1])
    ylim([0 3.5])
    set(gca,'xtick',[0:0.2:1],'xticklabels',100*[0:0.2:1],'ytick',[0:3],'fontsize',18,'ticklength',0.01*[1 1])

    set(fh_corr,'PaperPosition',[0 0 7 5],'PaperSize',[7 5]); 
end