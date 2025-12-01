function decoded_all= RMSCDecoding(r, m, codewords, p)
    [Nwords, N] = size(codewords);   % many codewords
    decoded_all= zeros(Nwords, N);           % result matrix
    mondegs = m - sum(de2bi(0:N-1, m), 2);
    A = find(mondegs <= r);     % 1-based info bit positions

    for i = 1:Nwords
        progressbar(i, Nwords)
        decoded_all(i, :) = RMSCDec(m,A,p, codewords(i, :));
    end
end
function u_hat=RMSCDec(m,A,p,codeword)
    N=2^m;            
    frozen_mask = true(N,1);
    frozen_mask(A) = false;
    y_probs = zeros(N,2);
    y_probs(:,1) = (1-p).*(codeword==0) + p.*(codeword==1);  
    y_probs(:,2) = p.*(codeword==0)   + (1-p).*(codeword==1);
    P = zeros(m + 1, N, 2); 
    B = zeros(m + 1, N);
    P(1, :, 1) = y_probs(:, 1); 
    P(1, :, 2) = y_probs(:, 2);
    for phi = 0 : (N - 1)
        recursivelyCalcP(m, phi);
        idx_m = phi + 1; 
        if frozen_mask(phi + 1)
            B(m + 1, idx_m) = 0;
        else
            prob_0 = P(m + 1, idx_m, 1);
            prob_1 = P(m + 1, idx_m, 2);
            
            if prob_0 > prob_1
                B(m + 1, idx_m) = 0;
            else
                B(m + 1, idx_m) = 1;
            end
        end
        if mod(phi, 2) == 1
            recursivelyUpdateB(m, phi);
        end
    end
    u_hat = B(m+1, :)'; 
    function recursivelyCalcP(lambda, phi)
        if lambda == 0
            return;
        end
        psi = floor(phi / 2);
        if mod(phi, 2) == 0
            recursivelyCalcP(lambda - 1, psi);
        end
        
        step = 2^(m - lambda);
        limit = step - 1; % 0 to 2^(m-lambda) - 1
        beta_range = 0:limit;
        idx_curr = (phi * step) + beta_range + 1;
   
        step_prev = 2^(m - (lambda - 1));
        idx_prev_0 = (psi * step_prev) + (2 * beta_range) + 1;     
        idx_prev_1 = (psi * step_prev) + (2 * beta_range + 1) + 1; 
        
        p_prev_0_0 = P(lambda, idx_prev_0, 1); % P(u' = 0)
        p_prev_0_1 = P(lambda, idx_prev_0, 2); % P(u' = 1)
        p_prev_1_0 = P(lambda, idx_prev_1, 1); % P(u'' = 0)
        p_prev_1_1 = P(lambda, idx_prev_1, 2); % P(u'' = 1)
        
        if mod(phi, 2) == 0
            P(lambda + 1, idx_curr, 1) = 0.5 * (p_prev_0_0 .* p_prev_1_0 + p_prev_0_1 .* p_prev_1_1);
            P(lambda + 1, idx_curr, 2) = 0.5 * (p_prev_0_1 .* p_prev_1_0 + p_prev_0_0 .* p_prev_1_1);
            
        else
            idx_neighbor = ((phi - 1) * step) + beta_range + 1;
            u_prime = B(lambda + 1, idx_neighbor); 
            p_upper_for_0 = zeros(size(u_prime));
            p_upper_for_0(u_prime == 0) = p_prev_0_0(u_prime == 0);
            p_upper_for_0(u_prime == 1) = p_prev_0_1(u_prime == 1);
            
            p_upper_for_1 = zeros(size(u_prime));
            p_upper_for_1(u_prime == 0) = p_prev_0_1(u_prime == 0);
            p_upper_for_1(u_prime == 1) = p_prev_0_0(u_prime == 1);
            
            P(lambda + 1, idx_curr, 1) = 0.5 * p_upper_for_0 .* p_prev_1_0;
            P(lambda + 1, idx_curr, 2) = 0.5 * p_upper_for_1 .* p_prev_1_1;
        end
    end

    function recursivelyUpdateB(lambda, phi)
        psi = floor(phi / 2);
        step = 2^(m - lambda);
        beta_range = 0 : (step - 1);
        idx_curr_even = ((phi - 1) * step) + beta_range + 1;
        idx_curr_odd  = (phi * step) + beta_range + 1;
        step_prev = 2^(m - (lambda - 1));
        idx_prev_even = (psi * step_prev) + (2 * beta_range) + 1;
        idx_prev_odd  = (psi * step_prev) + (2 * beta_range + 1) + 1;
        bit_left  = B(lambda + 1, idx_curr_even);
        bit_right = B(lambda + 1, idx_curr_odd);
        B(lambda, idx_prev_even) = xor(bit_left, bit_right);
        B(lambda, idx_prev_odd)  = bit_right;
        
        if mod(psi, 2) == 1
            recursivelyUpdateB(lambda - 1, psi);
        end
    end
end