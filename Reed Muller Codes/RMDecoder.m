function RMDecoder(r, m, decoder, p, file_name)
    n = 2^m;
    fid = fopen(file_name, 'r');
    save_names = ["_ml", "_scd", "_mld"];
    file_name_save = extractBefore(file_name, strlength(file_name) - 3);
    filesave = fopen(file_name_save+"_out"+ save_names(decoder) + ".txt", "w");
    allBits = fread(fid, 'ubit1');
    received_codewords = reshape(allBits, [size(allBits, 1)/n, n]);
    if decoder==1
        decoded_codewords = RMMLDecoderSimple(r, m, received_codewords);
    elseif decoder == 2
        decoded_codewords = RMSCDecoding(r, m, received_codewords, p);
    elseif decoder == 3
        decoded_codewords = RMMajorityDecoder(r, m, received_codewords);
    else
        disp("Invalid")
    end

    fwrite(filesave, decoded_codewords, 'ubit1');
    fclose(filesave);
end
