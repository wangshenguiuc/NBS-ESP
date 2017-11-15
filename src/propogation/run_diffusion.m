% [Input]
% A: adjacency matrix (could be weighted)
% method: 'personalized-pagerank' or 'self-diffusion'
% aux: auxiliary parameters, such as max number of iterations and
%      restart probability
%
% [Output]
% Q: diffusion state matrix. i-th column represents the diffusion
%    state of the i-th node.
%
function [Q] = run_diffusion(A, rst_p, reset,converge_rate)
n = size(A, 1);

renorm = @(M) bsxfun(@rdivide, M, sum(M));

A = A + diag(sum(A) == 0); % Add self-edges to isolated nodes
P = renorm(A)';
% reset =  renorm(reset_mat);
if ~exist('converge_rate','var')
    converge_rate = 1e-4;
end
norm_reset = renorm(reset')';
Q = norm_reset;
for i = 1:50
    Q_new = rst_p * norm_reset + (1 - rst_p) *Q*P;
    delta = norm(Q - Q_new, 'fro');
%     if mod(i,1)==0
%     fprintf('Iter %d. Frobenius norm: %f\n', i, delta);
%     end
    Q = Q_new;
    if delta < converge_rate
%         fprintf('Converged.\n');
        break
    end
end

% 
% Q = bsxfun(@rdivide, Q, sum(Q));

end
