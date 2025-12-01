
function PolarDecoder(m, A_file_name, decoder, p, file_name)
    n = 2^m;
    fid = fopen(file_name, 'r');
    save_names = ["_ml", "_scd"];
    file_name_save = extractBefore(file_name, strlength(file_name) - 3);
    filesave = fopen(file_name_save+"_out"+ save_names(decoder) + ".txt", "w");
    allA = matfile(A_file_name);
    A = allA.A;
    allBits = fread(fid, 'ubit1');
    received_codewords = reshape(allBits, [size(allBits, 1)/n, n]);
    if decoder==1
        decoded_codewords = PolarMLDecoder(m, A, received_codewords);
    elseif decoder==2
        decoded_codewords = PolarSCDecoder(m, A, p, received_codewords);
    end

    fwrite(filesave, decoded_codewords, 'ubit1');
    fclose(filesave);
    
end

function G = Generator(m)
     n = 2^m;
     G=zeros(n, n);
     for i=1:1:n
         ui = zeros(1,n);
         ui(i)=1;
         G(i,:)=RMmmencode(ui, m);
     end
end