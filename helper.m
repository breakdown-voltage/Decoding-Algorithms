 function G = Generator(m)
     n = 2^m;
     G=zeros(n, n);
     for i=1:1:n
         ui = zeros(1,n);
         ui(i)=1;
         G(i,:)=RMmmencode(ui, m);
     end
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