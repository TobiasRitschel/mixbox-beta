function Mat = legendre_to_bernstein(q)

% Extract parameters
M = q(1);

% Allocate memory
Mat = zeros(M+1, M+1);

for j = 0:M
    for k = 0:M
        % Initialize
        aux = 0;

        for i = max(0, j+k-M):min(j, k)
            aux = aux + (-1)^(k+i)*nchoosek(j, i)*nchoosek(k, i)*nchoosek(M-j, k-i);
        end

        % Compute element
        Mat(j+1, k+1) = aux/nchoosek(M, k);
    end
end