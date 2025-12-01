function recovered_codewords = RMMLDecoder(r, m, codewords)
    N = size(codewords, 1);
    n = 2^m;
    [G, A] = RMgenerator(r, m);
    G = gf(G(A, :), 1);
    [k, ~] = size(G);
    min_dist = n*ones(N);
    recovered_codewords = gf(zeros(N, n), 1);
    K = zeros(2^k, k);
    K(1, :) = dec2bin(0, k) - '0';
    for i=(1):1:((2^(k))-1)
        K(i+1, :) = dec2bin(i, k) - '0';
    end
    K = gf(K, 1);
    codespace = K * G;  %(2^x, n) codewords
    % size(codespace)
    for i=1:1:N
        codeword = codewords(i, :);
        codeword = repmat(codeword, 2^k, 1); %Transforming it to (2^k, k) to perform operations
        codeword = gf(codeword, 1);
        X = codespace + codeword; %Difference of the vectors 
        dist = sum(X~=0, 2); %Hamming distance of difference
        [val, argmin] = min(dist); %Minimum distance code
        if val<min_dist(i)
            recovered_codeword = codespace(argmin, :);
        end
        % disp(argmin)
        recovered_codewords(i, :) = recovered_codeword;
    end
    recovered_codewords = recovered_codewords.x;
end

function [G, A] = RMgenerator(r, m)
    n = 2^m;
    G=zeros(n, n);
    for i=1:1:n
        ui = zeros(1,n);
        ui(i)=1;
        G(i,:)=RMmmencode(ui, m);
    end
    mondegs = m-sum(de2bi(0:n-1, m), 2);
    A = find(mondegs <= r);
end