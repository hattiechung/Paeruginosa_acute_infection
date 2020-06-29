function [psamples, days_from_start, sampletypes] = get_patient_sample_dates(pat)

    sampleimport = importdata('../sample_dates.xlsx'); 
    samplenames = sampleimport.textdata(2:end,1); 
    sampledates = sampleimport.data(:,1);
    
    [~,pidxbinary] = get_patient_samples(pat,samplenames); 
    psamples = strtrim(samplenames(pidxbinary)); 
    pdates = sampledates(pidxbinary); 
    
    days_from_start = pdates - pdates(1) + 1;
    
    sampletypes = zeros(numel(psamples),1); 
    sorted_psamples = sort_nat(psamples);
    for s = 1:numel(sorted_psamples) 
        sampletypes(find(~cellfun(@isempty,strfind(psamples,sorted_psamples{s})))) = s; 
    end
    
    [days_from_start,sortorder] = sort(days_from_start,'ascend'); 
    psamples = psamples(sortorder); 
    sampletypes = sampletypes(sortorder); 
    
end