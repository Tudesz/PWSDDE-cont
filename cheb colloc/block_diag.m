function A_b = block_diag(A,square)
%BLOCK_DIAG Convert array of flattened out matricies to a Block diagonal
%matrix (for vector space and event function Jacobian matrix formulation)
% Input:
%   A: M x (n*n) or M x (1*n) matrix, each row containing a flattened n*n
%   or 1*n matrix
%   square: if true assume an n*n matrix otherwise 1*n (default true)
% Output:
%   A_b: (M*n) x (M*n) matrix containing Block diagonal entries of based on 
%       the columns of A

if nargin<2
    square = true;
end
% Initialization
M = size(A,1);
if square
    n = sqrt(size(A,2));
    A_b = zeros(M*n);
else
    n = size(A,2);
    A_b = zeros(M,M*n);
end

% Fill up matrix block by block with diagonal entries
for i = 1:n
    j_i = (i-1)*M+1:i*M; % row indicies (size M)
    if square
        for ii = 1:n
            j_ii = (ii-1)*M+1:ii*M; % column indicies (size M)
            A_b(j_i,j_ii) = diag(A(:,(ii-1)*n+i));
        end
    else
        A_b(:,j_i) = diag(A(:,i)); % only use column indicies (size M)
    end
end

end

