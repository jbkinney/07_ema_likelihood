function nseq = seq2nseq(seq)
    %nseq = double(1*(seq=='a' | seq=='A') + 2*(seq=='c' | seq=='C') + 3*(seq=='g' | seq=='G') + 4*(seq=='t' | seq=='T') + 5*(seq=='-'));
    nseq = double(1*(seq=='a' | seq=='A') + 2*(seq=='c' | seq=='C') + 3*(seq=='g' | seq=='G') + 4*(seq=='t' | seq=='T'));
end