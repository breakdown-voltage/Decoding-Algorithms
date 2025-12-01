function recovered_codewords = RMMajorityDecoder(r, m, codeword)
    n = 2^m;
    [G, ~] = RMgenerator(r, m);
    permutation_matrix = permutation(m);
    dec_matrix = idx_map_matrix(permutation_matrix);
    N = size(codeword, 1);
    recovered_codewords = gf(zeros(N, n), 1);
    for l=1:1:N
        recovered_codeword = gf(codeword(l, :), 1);
        %Take all the codewords with same degree
        for deg=r:-1:0
            mon_indices = sum(permutation_matrix, 2) == deg;
            mon_terms = permutation_matrix(mon_indices, :);
            num_terms = size(mon_terms, 1);
            check_dim = m - deg;
            check_vars_permutation = permutation(check_dim);
            for i=1:1:num_terms
                mon = mon_terms(i, :);
                check_vars = find(mon==0); %gets indices of the vars
                maj_arr = zeros(1, 2^check_dim);
                for j=1:2^(check_dim)
                    check_sum = gf(0, 1);
                    fixed_vals = check_vars_permutation(j, :);
                    [row_indices, num_indices] = get_indices(check_vars, fixed_vals, permutation_matrix);
                    for k=1:num_indices
                        index = dec_matrix(row_indices(k));
                        check_sum = check_sum + recovered_codeword(1, index);
                    end
                    maj_arr(j) = double(check_sum.x);
                end
                if sum(maj_arr == 1)>sum(maj_arr == 0)
                    % sprintf("1 Sum: %d| sum: %d", sum(maj_arr == 1), sum(maj_arr == 0))
                    target_mon_row = getidx(mon, m, G);
                    recovered_codeword = gf(recovered_codeword + target_mon_row, 1);
                end
            end
        end
        recovered_codewords(l, :) = recovered_codeword + gf(codeword(l, :), 1);
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

function row = getidx(a, m, G)
    row_idx = 1;
    for i=0:m-1
        if a(i+1)==0
            row_idx = row_idx + 2^i;
        end
    end
    row = gf(G(row_idx, :), 1);
end

function K = permutation(k)
    K = zeros(2^k, k);
    K(1, :) = dec2bin(0, k) - '0';
    for i=(1):1:((2^(k))-1)
        K(i+1, :) = dec2bin(i, k) - '0';
    end
end

function [row_indices, num_indices] = get_indices(cols, vals, matrix)
    temp = matrix(:, cols);
    operator = ismember(temp, vals, "rows");
    row_indices = find(operator);
    num_indices = size(row_indices, 1);
end
function x = idx_map_matrix(A)
    N = size(A, 1);
    x = zeros(N, 1);
    for i=1:1:N
        vector = A(i, :);
        x(i) = binaryVectorToDecimal(vector, 'MSBFirst') + 1;
    end
end
