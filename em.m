% FUNCTION NAME:
%   em
%
% DESCRIPTION:
%   Expectation maximization (EM) algorithm, which reconstructs a hyperspectral cube f from a computed tomgoraphy imaging spectrometer (CTIS) image g
%
% INPUT:
%   H       - System matrix H for the linear imaging equation g = Hf
%   g       - CTIS image
%   iter    - # of EM iterations
%
% OUTPUT:
%   f - Reconstructed hyperspectral cube
%
% ASSUMPTIONS AND LIMITATIONS:
%
%
% REVISION HISTORY
%   21/3/2022 - Mads Svanborg Peters
%       * Initial implementation
%

function f = em(H,g,iter)
g = g(:);           % Vectorisation of g
Ht = H.';           % Transposition of H, H^T
fk = Ht*g;          % Initial guess, alternatively, fk = ones(size(H,2),1);

hsum = sum(H).';    % Summation of rows in H
for i = 1:iter
    g_est = H*fk;                                   % Estimated CTIS image based on hyperspectral cube f
    g_ratio = g./g_est;                             % Ratio between true and estimated CTIS image
    g_ratio(isnan(g_ratio) | isinf(g_ratio)) = 0;   % Set inf and NaN values to zero
    fk_norm = fk./hsum;                             % Estimated cube fk normalized by hsum
    fk_norm(isnan(fk_norm) | isinf(fk_norm)) = 0;   % Set inf and NaN values to zero
    fk = fk_norm .* (Ht*g_ratio);                   % Update step to determine new guess fk+1
    disp(['EM iteration ',num2str(i),' out of ',num2str(iter)]);
end
f = fk;
end