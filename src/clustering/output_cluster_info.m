function output_cluster_info(indClust,nclst,pat_map_rev,pat_survival,sur_type,method_prefix)
valid_cluster = 0;
for i=1:nclst
    if length(find(indClust==i))>0
        valid_cluster = valid_cluster+1;
    end
    %     fprintf('cls %d:%d\t',i,length(find(indClust==i)));
end
% fprintf('\n')
fout = fopen(method_prefix,'w');
fprintf(fout,'name\ttime\tevent\tcluster\n');
for p = 1:length(indClust)
    pname = pat_map_rev(p);
    if isempty(pat_survival)
        fprintf(fout,'%s\t%f\t%d\t%d\n',char(pname),0,0,indClust(p));
        
    else
        if ~isKey(pat_survival.time{sur_type},pname)
%             fprintf('no key %d\n',p);
            continue;
        end
        if isnan(indClust(p))
%             fprintf('NAN %d\n',p);
            continue
        end
        fprintf(fout,'%s\t%f\t%d\t%d\n',char(pname),pat_survival.time{sur_type}(pname),pat_survival.event{sur_type}(pname),indClust(p));
    end
end
fclose(fout);

end