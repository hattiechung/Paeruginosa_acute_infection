patient_keys = {'05','13','14','09','45','42','01','46','10'};
patient_labels = {'A','B','C','D', 'E*', 'F*','G*','H*','I*'};
patient_conversion = containers.Map(patient_keys,patient_labels); 

patients_letter = {}; 
for pl = 1:numel(patient_labels)
    patients_letter{end+1} = patient_conversion(patient_keys{pl}); 
end

%%
data = importdata('mutation_NS_master_list_with_indels_20190607.xlsx'); 
genes = unique(data(2:end,4),'stable'); 
genes = genes(~cellfun(@isempty,genes));

% create gene-to-row index dict 
gene_idx = containers.Map(); 
for i = 1:numel(genes)
    gene_idx(genes{i}) = i;
end

% create patient-to-col index dict 
pat_idx = containers.Map(); 
for p = 1:numel(patients)
    pat_idx(patients{p}) = p; 
end
% initialize table

mut_cts = zeros(numel(genes),numel(patients));
for m = 2:size(data)
    % get gene 
    mut = data(m,3); 
    gene_entry = data(m,4); 
    gene = gene_entry{1};
    
    if ~strcmp(gene,'')
        % get row number for gene
        m_row = gene_idx(gene); 

        % get patient number for gene 
        ltag = data(m,2); 
        patient = ltag{1}(6:7); % THIS IS A HACK 
        m_col = pat_idx(patient); 

        % add to count
        mut_cts(m_row,m_col) = mut_cts(m_row,m_col)+1; 
    end
end

[sorted_cts,sort_idx]=sort(sum(mut_cts,2),'descend'); 
sort_idx = find(sum(mut_cts,2)>=2); 
mut_cts_mm = mut_cts(sort_idx,:);

% reorder gene names
mm_gene_names = genes(sort_idx);
mm_gene_names = strrep(mm_gene_names,'P01_unk','unk'); 
mm_ordered_names = fliplr({'pvdS','lasR','ladS','vfr','bifA','kinB',... %regulators
                    'mexR','ampD','sltB1',... %resistance
                    'unk'});
mm_idx = []; 
for mi = 1:numel(mm_ordered_names)
    for mj = 1:numel(mm_gene_names)
        if strcmp(mm_ordered_names{mi},mm_gene_names{mj})
            mm_idx = [mm_idx mj];
        end
    end
end

fh = figure; hold on; 
set(fh,'PaperSize',5*[1 1], 'PaperPosition',5*[0 0 1 1]); 
imagesc(mut_cts_mm(mm_idx,:)); 
set(gca,'ytick',1:numel(mm_gene_names(mm_idx)),'yticklabel',mm_gene_names(mm_idx), ...
    'xtick',1:numel(patients_letter),'xticklabel',patients_letter, ...
    'fontsize',18,'ticklength',0*[1 1]); 
% xlabel('Patients','fontsize',24)
% ylabel('Genes under selection','fontsize',24)

cmap = ([1 1 1; cbrewer('seq','Greys',4)]).^3;
colormap(cmap); 
axis equal; 
ylim([0.5 numel(sort_idx)+0.5])
xlim([0.5 numel(patients_letter)+0.5])

% manually plot grid
for x = 0:9
    plot([x+0.5,x+0.5],[0,10.5],'color',0.7*[1 1 1],'linewidth',0.3)
end
for y = 0:10
    line([0,10.5],[y+0.5,y+0.5],'color',0.7*[1 1 1],'linewidth',0.3)
end
figdir = '~/Dropbox/projects/kishony-acute-infection/analysis_assembly/figures'; 
hold off; 
saveas(fh,[figdir filesep 'genes_under_selection_table_indels.png'],'png');

fb_cb = figure; 
set(fb_cb,'PaperSize',6*[1 1],'PaperPosition',6*[0 0 1 1]); 
cb = colorbar;
colormap(cmap); 
axis off; 
lim = 0:4; 
caxis([min(lim) max(lim)])
set(cb, 'Ytick', lim, 'yticklabel', lim)
set(cb,'position',[.15 .1 .04 .3]);
set(cb, 'fontsize', 20)
% saveas(fb_cb,[figdir filesep 'genes_under_selection_colorbar_indels.png'],'png');

cgo = clustergram(mut_cts,'ColumnLabels',patients_letter, 'RowLabels',genes, 'Colormap', 'redbluecmap');