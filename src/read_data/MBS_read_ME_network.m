function [ all_pv_net,pat_diff] = MBS_read_ME_network( cancer_type,  rst, p2m, i2hgnc, topk, ME_hop_l)
%READ_ME_NETWORK Summary of this function goes here
%   Detailed explanation goes here
[npat,ngene] = size(p2m);
mute_freq = sum(p2m);
if ~exist('ME_hop_l','var')
    ME_hop_l = 1:3;
end
pv_thres_init = 0.5;
pv_thres = pv_thres_init;
for ME_hop = ME_hop_l
    ME_file =  ['../data/ME_network/',cancer_type,'_ME_',num2str(ME_hop),'.txt'];
    if ~exist(ME_file,'file')
        continue
    end
     [gid1,gname1,gid2,gname2,pv,odds] = textread(ME_file,'%d%s%d%s%f%f');
    nnet_edge =length(gid1);
    pv_left_mat = ones(ngene,ngene);
    for ii=1:nnet_edge
        pv_left_mat(gid1(ii),gid2(ii)) = pv(ii);
    end    
    [vi,vj] = find(pv_left_mat<pv_thres);
    nedge = length(vi);
    if nedge>=100
        break
    end
end

[o,v] = sort(pv_left_mat(:),'ascend');
topk = min(topk,length(o));
pv_thres = min(pv_thres_init,o(topk));
all_pv_net = pv_left_mat <= pv_thres;
% ME_network_dir = '../inter_result/ME_network/';
% if ~exist(ME_network_dir,'dir')
%     mkdir(ME_network_dir)
% end
% fout = fopen([ME_network_dir,cancer_type,num2str(pv_thres),'_',num2str(topk),'.txt'],'w');
% [vi,vj] = find(all_pv_net);
% for i=1:length(vi)
%     fprintf(fout,'%s\t%s\t%f\n',i2hgnc(vi(i)),i2hgnc(vj(i)),pv_left_mat(vi(i),vj(i)));
% end
% fclose(fout);
pat_diff = run_diffusion(all_pv_net,rst,p2m);
nan_pat = isnan(sum(pat_diff,2));
pat_diff(nan_pat,:) = 0;
end

