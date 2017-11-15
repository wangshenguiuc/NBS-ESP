function [pat_mute_mat,pat_survival,network,pat_map_rev,gene_map_rev,pat_map,gene_map] = ...
    MBS_read_data(mutation_file,network_file,survival_file,min_mutation,weighted,filter_no_survival_pat,read_CNA)
%READ_DATA Summary of this function goes here
%   Detailed explanation goes here

if ~exist('read_CNA', 'var')
    read_CNA = false;
end

[p,t1,e1,t2,e2,stg] = textread(survival_file,'%s%f%d%f%d%s','emptyvalue',-1,'headerlines',1,'delimiter','\t');
filter = (t1~=-1)&(t2~=-1);
p = p(filter);
survial_p = unique(p);

if ~read_CNA
[pat,mut_gene] = textread(mutation_file,'%s%s');
else
    [pat,mut_gene,isCNA] = textread(mutation_file,'%s%s%s');
    filter = strcmp(isCNA,'CNA');
    pat = pat(filter);
    mut_gene = mut_gene(filter);
end
if filter_no_survival_pat
    pat_set = intersect(unique(pat),survial_p);
else
    pat_set = unique(pat);
end
if isempty(pat_set) || length(pat_set)==1
    pat_mute_mat=0;
    pat_survival=[];
    network=0;
    pat_map_rev=0;
    gene_map_rev=0;
    gene_map = 0;
    pat_map=0;
    return
end



gene_set = unique(mut_gene);

ngene = length(gene_set);
gene_map = containers.Map(gene_set,1:ngene);
gene_map_rev = containers.Map(1:ngene,gene_set);
npat = length(pat_set);
tmp_pat_map = containers.Map(pat_set,1:npat);
tmp_pat_map_rev = containers.Map(1:npat,pat_set);

filter = isKey(tmp_pat_map,pat);
pat = pat(filter);
mut_gene = mut_gene(filter);

pid = cell2mat(values(tmp_pat_map,pat));
gid = cell2mat(values(gene_map,mut_gene));
tmp_pat_mute_mat = sparse(pid,gid,1,npat,ngene);

pat_valid = find(sum(tmp_pat_mute_mat,2)>=min_mutation);
pat_valid_id = values(tmp_pat_map_rev,num2cell(pat_valid));
pat_set = unique(pat_valid_id);
npat = length(pat_set);
pat_map = containers.Map(pat_set,1:npat);
pat_map_rev = containers.Map(1:npat,pat_set);

if ~read_CNA
[pat,mut_gene] = textread(mutation_file,'%s%s');
else
    [pat,mut_gene,isCNA] = textread(mutation_file,'%s%s%s');
    filter = strcmp(isCNA,'CNA');
    pat = pat(filter);
    mut_gene = mut_gene(filter);
end

filter = isKey(pat_map,pat);
pat = pat(filter);
mut_gene = mut_gene(filter);
gene_set = unique(mut_gene);
ngene = length(gene_set);
gene_map = containers.Map(gene_set,1:ngene);
gene_map_rev = containers.Map(1:ngene,gene_set);

pid = cell2mat(values(pat_map,pat));
gid = cell2mat(values(gene_map,mut_gene));
pat_mute_mat = full(sparse(pid,gid,1,npat,ngene));

if weighted
    [e1,e2,r] = textread(network_file,'%s%s%f','delimiter','\t');
else
    [e1,e2] = textread(network_file,'%s%s','delimiter','\t');
    r = ones(size(e1));
end
filter = isKey(gene_map,e1) & isKey(gene_map,e2);
r = r(filter);
e1 = e1(filter);
e2 = e2(filter);
e1 = cell2mat(values(gene_map,e1));
e2 = cell2mat(values(gene_map,e2));
filter = e1<e2;
e1 = e1(filter);
e2 = e2(filter);
r = r(filter);
network = sparse(e1,e2,r,ngene,ngene);
network(network>1) = 1;
network = max(network,network');

[p,t1,e1,t2,e2,stg] = textread(survival_file,'%s%f%d%f%d%s','emptyvalue',-1,'headerlines',1,'delimiter','\t');
filter = (t1~=-1)&(t2~=-1);
p = p(filter);
t1 = t1(filter);
e1 = e1(filter);
t2 = t2(filter);
e2 = e2(filter);
if length(t1)==0
    pat_survival = [];
else
    filter = isKey(pat_map,p);
    p = p(filter);
    t1 = t1(filter);
    e1 = e1(filter);
    t2 = t2(filter);
    e2 = e2(filter);
    pat_survival.time{1} = containers.Map(p,t1);
    pat_survival.event{1} = containers.Map(p,e1);
    pat_survival.time{2} = containers.Map(p,t2);
    pat_survival.event{2} = containers.Map(p,e2);
end

end

