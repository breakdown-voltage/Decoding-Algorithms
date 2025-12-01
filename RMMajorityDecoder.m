function recovered_codewords = RMMajorityDecoder(r, m, codeword)
    n = 2^m;
    [G, ~] = RMgenerator(r, m);
    permutation_matrix = permutation(m);
    dec_matrix = idx_map_matrix(permutation_matrix);
    N = size(codeword, 1);
    recovered_codewords = gf(zeros(N, n), 1);
    for l=1:1:N
        recovered_codeword = gf(codeword(l, :), 1);
        progressbar(l, N);
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