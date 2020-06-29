function [patientidx,patientidxbinary] = get_patient_samples(patient,SampleNames)
    matches = strfind(SampleNames,patient); 
    patientidx = []; 
    patientidxbinary = logical(zeros(1,numel(SampleNames))); 
    for m = 1:numel(matches) 
        cur_match = matches{m}; 
        if ~isempty(cur_match)
            if cur_match==1
                patientidx = [patientidx m]; 
                patientidxbinary(m) = 1; 
            end
        end
    end
end