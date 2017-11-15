function [ g2hgnc,i2hgnc,hgnc2i,hgnc2g] = map_gene_name(g2i,i2g,gene_name_mapping_file)
%MAP_GENE_NAME Summary of this function goes here
%   Detailed explanation goes here
[entrez,hgnc] = textread(gene_name_mapping_file,'%s%s','headerlines',1);

ngene = length(g2i);

gname = values(i2g,num2cell(1:ngene));

[Lia,Lib] = ismember(entrez,gname);

entrez = entrez(Lia);
hgnc = hgnc(Lia);
gid = Lib(Lia);

g2hgnc = containers.Map(entrez,hgnc);
hgnc2g = containers.Map(hgnc,entrez);
i2hgnc = containers.Map(gid,hgnc);
hgnc2i = containers.Map(hgnc,gid);


end

