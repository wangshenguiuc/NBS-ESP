function [ res_pv,min_nclst,max_nclst] = get_survival_pvalue( pat_survival,file_name,sur_type,indClust,pat_map_rev,nclst,write_to_file)
%GET_SURVIVAL_PVALUE Summary of this function goes here
%   Detailed explanation goes here
if sur_type==1
    sur_type_str = 'ov';
else
    sur_type_str = 'dfs';
end
indClust = full(indClust);
method_prefix = [ file_name sur_type_str '.pvalue'];
cls_index_file = [file_name sur_type_str '.cluster'];
fout = fopen(cls_index_file,'w');
remove_c = [];
for i=1:nclst
    cct = length(find(indClust==i));
    if cct<1
        remove_c = [remove_c,i];
    end
end
max_nclst = -1;
min_nclst = 10000;
for i=1:nclst
    if length(find(indClust==i)) > max_nclst
        max_nclst = length(find(indClust==i));
    end
    if length(find(indClust==i)) < min_nclst
        min_nclst = length(find(indClust==i));
    end
end
% fprintf('clst:%d %d %d %d\n',nclst,i,max_nclst,min_nclst);

fprintf(fout,'time\tevent\tour_cluster\tnbs_cluster\n');
time = [];
clst = [];
event = [];
clst_list = [];
for p = 1:length(indClust)
    pname = pat_map_rev(p);
    if ~isKey(pat_survival.time{sur_type},pname)
        continue;
    end
    if isnan(indClust(p))
        continue
    end
    if ~isempty(find(remove_c==indClust(p), 1))
        continue
    end
    time = [time,pat_survival.time{sur_type}(pname)];
    event = [event,pat_survival.event{sur_type}(pname)];
    clst = [clst,indClust(p)];
    clst_list = [clst_list,indClust(p)];
    fprintf(fout,'%f\t%d\t%d\t%d\n',pat_survival.time{sur_type}(pname),pat_survival.event{sur_type}(pname),indClust(p),indClust(p));
end
fclose(fout);
min_csize = min(length(find(clst==1)),length(find(clst==2)));
% if max(clst)==2 && min_csize<=5
%     x1 = [time(clst==1),event(clst==1)];
%     x2 = [time(clst==2),event(clst==2)];
%     [chisq,pv] = logrank(x2,x1);
%     res_pv = inf(1,12);
%     res_pv(8) = chisq;
%     res_pv(4) = pv;
%     return
% end

clst_list = unique(clst_list);


ts_l = [60,120,150,180];

[x, y] =system(['Rscript --vanilla evaluating_survival/NBS_survival_prediction.R ', ...
    cls_index_file,' ',method_prefix,' 1',' ',num2str(ts_l(1)),' ',...
    num2str(ts_l(2)),' ',num2str(ts_l(3)),' ',num2str(ts_l(4)),' ',num2str(nclst)-length(remove_c)]);

if length(clst_list)==1
    res_pv = ones(1,12);
    return
end
if x==1
    error('r script error %s',y);
    
else
    res_pv = dlmread(method_prefix);
end

if exist('write_to_file','var')
    
    fout = fopen(cls_index_file,'w');
    
    fprintf(fout,'time\tevent\tcluster\tnbs_cluster\n');
    time = [];
    clst = [];
    event = [];
    for p = 1:length(indClust)
        pname = pat_map_rev(p);
        if ~isKey(pat_survival.time{sur_type},pname)
            continue;
        end
        if isnan(indClust(p))
            continue
        end
        if ~isempty(find(remove_c==indClust(p), 1))
            continue
        end
        time = [time,pat_survival.time{sur_type}(pname)];
        event = [event,pat_survival.event{sur_type}(pname)];
        clst = [clst,indClust(p)];
        
        fprintf(fout,'%f\t%d\t%d\t%d\n',pat_survival.time{sur_type}(pname),pat_survival.event{sur_type}(pname),indClust(p),indClust(p));
    end
    fclose(fout);
    [x, y] =system(['Rscript --vanilla evaluating_survival/plot_our_km_plot.R ', ...
        cls_index_file,' ',write_to_file,' ',num2str(nclst)]);
    
    if x~=0
        error('r script error %s',y);
    end
end
end

