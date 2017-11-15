function [clust_pv,indClust_all,chisq,cluster_max,cluster_min] = our_NBS_cluster(pat_mute,integrate_mat,pat_survival,rst_prob,selected_pat,i2p,degree,cancer_type,max_nclst,min_nclst,...
    low_d_type_l,input_mat_type_l,clst_type_l,dist_type_l,dim_l,method_name)
addpath 'evaluating_survival'
nmf_options.iter = 100;
nmf_options.gamma = 200;
nmf_options.tof = 1e-4;
nmf_options.dis = false;
nmf_options.distance = 'nnls';

ngene = size(integrate_mat{1},2);
% assert(ngene>1)
diff_mat_type = 3;
if ~exist('method_name', 'var')
    method_name = '';
end
if ~exist('min_nclst', 'var')
    min_nclst = 2;
end
if ~exist('max_nclst', 'var')
    max_nclst = 8;
end

input_mat_type = input_mat_type_l(end);

clst_type = clst_type_l(end);
cluster_max = zeros(1,max_nclst);
cluster_min = zeros(1,max_nclst);
clust_pv = ones(1,max_nclst);
chisq = zeros(1,max_nclst);
dist_type = dist_type_l(end);

dim = dim_l(end);
var_remove_l = [0];

min_sumD = inf(1,8);
consen_indClust =cell(1,8);
for i=1:8
    consen_indClust{i} = [];
end
npara = length(dim_l)*length(dist_type_l)*length(input_mat_type_l);
nmat = length(integrate_mat);
dist_mat = cell(nmat,npara);


for mati =1:nmat
    diff_mat = integrate_mat{mati};
    
    select_pat_mute = pat_mute(selected_pat,:);
    select_diff_mat = diff_mat(selected_pat,:);
    % neg_pat_info = pat_info(selected_pat,:);
    nneg_pat = length(selected_pat);
    select_i2p = containers.Map(1:nneg_pat,values(i2p,num2cell(selected_pat)));
    std_prop = std(select_diff_mat);
    
    npat = size(select_diff_mat,1);
    %     fprintf('mati=%d,',mati);
    ngene = size(select_diff_mat,2);
    valid_gene_set = 1:ngene;
    if ngene==1
        for nclst=min_nclst:max_nclst
            rng(1);
            indClust = kmeans(full(reshape(select_diff_mat,npat,1)),nclst,'replicates',100);
            consen_indClust{nclst} = [consen_indClust{nclst},indClust];
        end
        continue
    end
    ct = 1;
    
    for var_remove = var_remove_l
        if var_remove == 1
            std_cutoff = 0.2;
            [o,v] = sort(std_prop,'ascend');
            thres = o(floor(ngene*std_cutoff));
            high_std_node = find(std_prop>thres);
            valid_gene_set = intersect(find(degree>0),find(std_prop>thres));
        end
        for input_mat_type = input_mat_type_l
            if input_mat_type == 1
                mat = select_diff_mat(:,valid_gene_set);
                %              sdiff = repmat(sum(mat,2),1,size(mat,2));
                %             mat = mat./sdiff;
            elseif input_mat_type == 2
                mat = select_diff_mat(:,valid_gene_set);
                sdiff = repmat(sum(mat,2),1,size(mat,2));
                mat = mat./sdiff;
            elseif input_mat_type == 3
                mat = log(select_diff_mat(:,valid_gene_set)+1/(ngene*ngene))-log(1/(ngene*ngene));
            elseif input_mat_type == 4
                mat = log(select_diff_mat(:,valid_gene_set)+1/(ngene*ngene))-log(1/(ngene*ngene));
                sdiff = repmat(sum(mat,2),1,size(mat,2));
                mat = mat./sdiff;
            else
                error('unkonwn input mat type\n');
            end
            %         dlmwrite(['../inter_result/propogate_matrix/',cancer_type,'_',num2str(norm_type),'_',num2str(input_mat_type),'_',num2str(dim_i),'_',num2str(dist_type),'.txt'],mat,'delimiter','\t');
            for dim = dim_l
                if dim == 0
                    dim_i = floor(npat/2);
                else
                    dim_i = dim;
                end
                if dim == -1
                    [~,~,~,~,explained] = pca(full(mat));
                    for ii=1:10
                        fprintf('%d\t%f\n',ii,explained(ii));
                    end
                    dim_i = find(cumsum(explained)>30, 1 );
                    %                     plot(1:length(explained),cumsum(explained));
                    %                   fprintf('dim=%d\n',dim_i)
                end
                for low_d_type = low_d_type_l
                    if low_d_type == 1
                        [pat_U,pat_S,~] = svds(mat,dim_i);
                        low_d_mat = pat_U*pat_S.^0.5;
                    elseif low_d_type == 2
                        [A,Y]=nmfnnls(mat,floor(npat/2),nmf_options);
                        low_d_mat = A;
                    elseif low_d_type == 3
                        [pat_U,pat_S,~] = svds(mat,dim_i);
                        low_d_mat = pat_U*pat_S.^0.5;
                    elseif low_d_type == 4
                        [A,Y]=nmfnnls(mat,dim_i,nmf_options);
                        low_d_mat = A;
                    elseif low_d_type == 5
                        low_d_mat = mat;
                    else
                        error('unkonwn input mat type\n');
                    end
                    for dist_type = dist_type_l
                        if dist_type == 1
                            dist = squareform(pdist(low_d_mat,'cosine'));
                        elseif dist_type == 2
                            dist = squareform(pdist(low_d_mat));
                        elseif dist_type == 3
                            dist = squareform(pdist(low_d_mat,'cosine'));
                            dist = squareform(pdist(dist,'cosine'));
                        elseif dist_type == 4
                            dist = mat;
                        elseif dist_type == 5
                            dist = squareform(pdist(low_d_mat,'cosine'));
                            dist = squareform(pdist(dist,'cosine'));
                            dist = squareform(pdist(dist,'cosine'));
                        elseif dist_type == 6
                            dist = squareform(pdist(low_d_mat));
                            dist = squareform(pdist(dist,'cosine'));
                        elseif dist_type == 7
                            dist = squareform(pdist(low_d_mat,'cosine'));
                            dist = squareform(pdist(dist));
                        end
                        %                         dist = zscore(dist,[],2);
                        dist_mat{mati,ct} = dist;
                        ct = ct + 1;
                        
                    end
                end
            end
        end
    end
end
%
if nmat>1 || size(integrate_mat{1},2)~=1
    
    for pi = 1:npara
        dist = 0;
        for mati=1:nmat
            diff_mat = integrate_mat{mati};
            ngene = size(diff_mat,2);
            if ngene==1
                continue
            end
            dist = dist_mat{mati,pi};
            for clst_type = clst_type_l
                for nclst=min_nclst:max_nclst
                    if clst_type == 1
                        rng(1);
                        [indClust,~,sumD] = kmeans(full(dist),nclst,'replicates',100);
                    elseif clst_type == 2
                        rng(1);
                        [indClust,~,sumD] = kmeans(full(dist),nclst,'distance','cosine','replicates',100);
                    end
                    consen_indClust{nclst} = [consen_indClust{nclst},indClust];
                end
            end
            
        end
    end
end

%
% if nmat>1 || size(integrate_mat{1},2)~=1
%
%     for pi = 1:npara
%         dist = 0;
%         for mati=1:nmat
%             diff_mat = integrate_mat{mati};
%             ngene = size(diff_mat,2);
%             if ngene==1
%                 continue
%             end
%             if nmat > 1
%                 dist = dist + zscore(dist_mat{mati,pi},[],1);
%             else
%                 dist = dist_mat{mati,pi};
%             end
%         end
%         for clst_type = clst_type_l
%             for nclst=min_nclst:max_nclst
%                 if clst_type == 1
%                     rng(1);
%                     [indClust,~,sumD] = kmeans(full(dist),nclst,'replicates',100);
%                 elseif clst_type == 2
%                     rng(1);
%                     [indClust,~,sumD] = kmeans(full(dist),nclst,'distance','cosine','replicates',100);
%                 end
%                 consen_indClust{nclst} = [consen_indClust{nclst},indClust];
%             end
%         end
%     end
% end
subtype_folder = ['../output/subtype/',cancer_type,'/'];
if exist(subtype_folder)~=7
    mkdir(subtype_folder)
end
indClust_all = zeros(npat,max_nclst);

km_plot_folder = ['../output/km_plot/',cancer_type,'/'];
if exist(km_plot_folder)~=7
    mkdir(km_plot_folder)
end
%   if nmat==1 && size(integrate_mat{1},2)==1
%         dim = dim_l(1);
%         NBS_prefix = [cancer_type,'_',num2str(nclst),'_',num2str(rst_prob),'_',num2str(diff_mat_type),'_', ...
%             num2str(integrate_type),'_mutation_sum'];
%     else

NBS_prefix = [cancer_type,'_',num2str(nclst),'_',num2str(rst_prob),'_',num2str(diff_mat_type),'_',num2str(input_mat_type), ...
    '_',num2str(clst_type), ...
    '_',num2str(dim),'_',num2str(dist_type),method_name,'.txt'];
%     end
cluster_pv_folder = ['../output/NBS_comparision/',cancer_type,'/'];
if exist(cluster_pv_folder)~=7
    mkdir(cluster_pv_folder)
end
fout = fopen([cluster_pv_folder,NBS_prefix],'w');
for nclst=min_nclst:max_nclst
    
    
    indClust_prefix = [cancer_type;
    
    
    indClust = aggregate_cluster(consen_indClust{nclst});
    indClust_all(:,nclst) = indClust;
    
    cls_size = zeros(1,nclst);
    for nc=1:nclst
        cls_size(nc) = length(find(indClust==nc));
    end
    cluster_max(nclst) = max(cls_size);
    cluster_min(nclst) = min(cls_size);
    method_prefix = [subtype_folder, indClust_prefix];
    output_cluster_info(indClust,nclst,select_i2p,pat_survival,1,method_prefix);
    
    
    file_name = [km_plot_folder,indClust_prefix];
    if isempty(pat_survival)
        continue
    end
    [ res_pv_dfs,~,~] = get_survival_pvalue(pat_survival,[file_name,'_dfs'],2,indClust,select_i2p,nclst);
    [ res_pv_ov,~,~] = get_survival_pvalue(pat_survival,[file_name,'_ov'],1,indClust,select_i2p,nclst);
    clust_pv(nclst) = res_pv_ov(4);
    chisq(nclst) = res_pv_ov(8);
    
    fprintf(fout,'%s\t%d\t%d\t%f\t',cancer_type,dim,nclst,rst_prob);
    for r = 1:length(res_pv_dfs)/2
        fprintf(fout,'%e\t',res_pv_dfs(r));
    end
    for r = 1:length(res_pv_ov)/2
        fprintf(fout,'%e\t',res_pv_ov(r));
    end
    fprintf(fout,'\n');
end
fclose(fout);

