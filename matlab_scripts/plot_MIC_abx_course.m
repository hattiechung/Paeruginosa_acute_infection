function [patient_med_dates, days_fr_sample1, pat_cult_fr_sample1] = plot_MIC_abx_course(pat,makeplot)
    % INPUT pat is a numeric
    
    if nargin<2
        makeplot = true;
    end
    spreadsheet_drug_names = {'azithro','cipro','tobramycin','mero','cefep','ceftaz','ceftriax','piptaz','vanc','mero'};
    drug_rename = containers.Map(spreadsheet_drug_names,...
                {'Azithromycin','Ciprofloxacin','Tobramycin','Meropenem','Cefepime','Ceftazidime','Ceftriaxone','PipTaz','Vancomycin','Meropenem'});
            
    relevant_drugs = {'cefep','ceftaz','ceftriax','mero','piptaz','cipro','azithro'};
    
    rxfile = '../all_rx.xlsx';
    data = importdata(rxfile); 
    patient = data.data(:,1); 
    rx_date = data.data(:,3); 
    rx = data.textdata(2:end,4); 
    
    % get med dates for patient
    patient_meds = unique(rx(patient==pat));
    unique_dates = unique(rx_date(patient==pat)); 
    patient_dates = min(unique_dates):max(unique_dates);
    patient_med_dates = zeros(numel(relevant_drugs),numel(patient_dates)); 
    for m = 1:numel(relevant_drugs)
        if find(ismember(patient_meds,relevant_drugs{m}))
            med_dates = unique(rx_date(patient==pat & ismember(rx,relevant_drugs{m}))); 
            patient_med_dates(m,ismember(patient_dates,med_dates))=1;
        end
    end
    
    % get sample dates for patient
    sampledata = importdata('../sample_dates.xlsx'); 
    culturenames = sampledata.textdata(2:end,1); 
    culturedates = sampledata.data(:,1);
    culturepatients = [];
    for c = 1:numel(culturenames)
        culturepatients = [culturepatients str2num(culturenames{c}(1:2))]; 
    end
    
    patcultdates = culturedates(culturepatients==pat);
    days_fr_sample1 = patient_dates - min(patcultdates);
    pat_cult_fr_sample1 = patcultdates - min(patcultdates);
    
    drugs_renamed = {};
    for d = 1:numel(relevant_drugs) 
        drugs_renamed = [drugs_renamed drug_rename(relevant_drugs{d})]; 
    end
    
    if makeplot
        fh = figure; hold on; 
        imagesc(patient_med_dates); 
        box off; 
        set(gca,'ytick',1:numel(relevant_drugs),'yticklabel',drugs_renamed,...
                'xtick',1:numel(days_fr_sample1),'xticklabel',days_fr_sample1,...
                'fontsize',12,'ticklength',0*[1 1]); 

        x_min = find(days_fr_sample1==0)-5.5; 
        x_max = length(days_fr_sample1)+0.5; 
        y_min = 0.5; 
        y_max = size(patient_med_dates,1)+0.5; 
        for x = x_min:x_max
            plot([x,x],[y_min-0.5,y_max+0.5],'color',0.8*[1 1 1],'linewidth',0.3)
        end
        for y = y_min-1:y_max+1
            line([x_min,x_max],[y,y],'color',0.8*[1 1 1],'linewidth',0.3)
        end
        hold off; 

        xtickangle(90)
        colormap(flipud(bone));
        axis equal;

        xlim([x_min x_max]);
        ylim([y_min y_max])

        set(fh,'PaperPosition',[0 0 6 2],'PaperSize',[6 2]);
    end
end