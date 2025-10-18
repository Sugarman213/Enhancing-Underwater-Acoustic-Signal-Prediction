function VSLDE = VS_Link_Dispersion_Entropy(x, m, c, t,n)
    % x: input time sequence vecto
    % m: embedding dimension
    % c: number of classes
    % d: time delay
    
    N = length(x);  % Length of the time sequence
    sigma_x = std(x);  % Standard deviation of x
    mu_x = mean(x);  % Mean of x
    
    % Step 1: Map the original sequence to a symbolic sequence by NCDF algorithm
%     y = normcdf(x, mu_x, sigma_x);  % NCDF mapping
%     y = round(y * (c - 1) + 0.5);  % Linear rounding to the nearest integer between 1 and c
%     z = mod(y, c) + 1;  % Ensure symbols are in the range [1, c]

%%   NCDF--Mapping approaches
    y=normcdf(x,mu_x,sigma_x);
    y=mapminmax(y,0,1);
    y(y==1)=1-1e-10;
    y(y==0)=1e-10;
    z=round(y*c+0.5);


%%

    % Step 2: Reconstruct the time sequences into a series of orbits
    L = N - (m - 1) * t;  % Adjusted length for phase-space reconstruction
    orbits = zeros(L, m);  % Initialize orbits matrix

    for i = 1:L
        for j = 1:m
            orbits(i, j) = z(i + (j - 1) * t);  % Reconstruct phase-space orbits
        end
    end
    
    % Step 3: Calculate the transition probability
    P = zeros(c^m);  % Initialize transition probability matrix
    for i = 1:L
        current_orbit = orbits(i, :);
        if i < L-n
            next_orbit = orbits(i + n, :);
            current_index = 1 + sum((current_orbit - 1) .* (c.^(0:(m - 1))));
            next_index = 1 + sum((next_orbit - 1) .* (c.^(0:(m - 1))));
            P(current_index, next_index) = P(current_index, next_index) + 1;
        end
    end
    % Normalize the rows of P to sum to 1
    P = P ./ repmat(sum(P, 2), 1, c^m);
    
    % Step 4: Calculate the probability of each type of link
    L_vect = zeros(1, c^m);  % Initialize L vectors
    for i = 1:c^m
        P_row = P(i, :);
        L_vect(i) = -sum(P_row(P_row > 0) .* log(P_row(P_row > 0)));  % Calculate entropy for each row
    end
    
    % Step 5: Calculate the LDE
    VSLDE = mean(L_vect);  % Calculate the LDE as the mean of L_vect 
end