function Lambda = beta_to_legendre(q)

% Extract parameters
M  = q(1);
dt = q(2);

% Allocate memory
Lambda = zeros(M+1, M+1);

for j = 0:M
    for k = 0:M
        % Initialize
        aux = 0;

        for i = 0:j
            aux = aux + (-1)^(j+i)*nchoosek(j, i)^2/nchoosek(M+j, k+i);
        end

        % Compute element
        Lambda(j+1, k+1) = (2*j + 1)/(M+j+1)*nchoosek(M, k)*aux;
    end
end

% Scale to interval [0, dt]
Lambda = Lambda*(M+1)/dt;