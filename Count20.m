%Return aligned glycopeptide sequences ranging from -20 to 20
%Uniprot entry for peptide sequence .xlsx is at column 2
%Glycosylation site within protein for peptide sequence .xlsx is at column 3
%Glycosylation site within peptide for peptide sequence .xlsx is at column
%4
%Uniprot entry for database .xlsx is at column 5
clear;
tb = readtable("Peptide.xlsx");
ptb = readtable("GlycoproteinFullLength_hilic1.xlsx");
%a: amino acid count (vector of 20)
%y: filename
%tb: table
%i1: number of row
%i3: number of column
%i4: count within vector b (used to identify aa)
%i5: count within vector for writing into output table
%i6: count within output table for number of row
%i7: count within output table for number of column
%i8: number of row in trembl table
%i9: number of glyposition retrieved from trembl
%i10: number of disagreement between pglyco and trembl
%dsg: row of disagreement between pglyco and trembl
outputtb = cell2table(cell(100000,41));
b = ['R';'H';'K';'D';'E';'S';'T';'N';'Q';'C';'G';'P';'A';'I';'L';'M';'F';'W';'Y';'V';'J'];
i1 = 1;
i4 = 1;
i5 = 1;
j = b(8);
p = b(12);
i6 = 1;
i7 = 1;
i8 = 1;
i9 = 0;
i10 = 0;
dsg = zeros(height(tb),1);
while i1 < (height(tb) + 1)
    i3 = table2array(tb(i1,4)) + 6;
    t = tb(i1,i3);
    r = table2array(t);
    e = cell2mat(r);
    if e == j
        t2 = tb(i1,(i3+1));
        r2 = table2array(t2);
        e2 = cell2mat(r2);
        if e2 ~= p
            i5 = i3 - 20;
            if i5 > 6
                %if there're at least 10 amino acids to the left of J
                i7 = 1;
                %initialize i7 after last round of counting
                while i7 < 42 && i5 < (width(tb) + 1)
                %write into output table
                    t = tb(i1,i5);
                    r = table2array(t);
                    e = cell2mat(r);
                    i4 = 1;
                    while i4 < 21
                        if e == b(i4)
                            outputtb(i6,i7) = tb(i1,i5);
                            i4 = 23;
                            %write into output table
                        end
                        i4 = i4 + 1;
                    end
                    i5 = i5 + 1;
                    i7 = i7 + 1;
                end
                i4 = 1;
                t = outputtb(i6,41);
                r = table2array(t);
                e = cell2mat(r);
                while i4 < 21
                    if e == b(i4)
                        i4 = 23;
                        %when the last aa in output table is filled up with
                        %aa, move to next row in output table, or stay in
                        %the same row and cover original data
                    end
                    i4 = i4 + 1;
                end
                if i4 == 21
                    entry1 = cell2mat(table2array(tb(i1,2)));
                    %entry of peptide whose J is to the end of sequence
                    prosite1 = tb(i1,3);
                    prosite3 = table2array(prosite1);
                    %glycosylation position
                    i8 = 1;
                    while i8 < (height(ptb) + 1)
                        entry2 = cell2mat(table2array(ptb(i8,5)));
                        if strcmpi(entry1,entry2)
                            gsite1 = ptb(i8,(prosite3 + 6));
                            gsite2 = table2array(gsite1);
                            gsite3 = cell2mat(gsite2);
                            if gsite3 == 'N'
                                i5 = prosite3 -14;
                                %position of the first aa from ptb to 
                                %be written into outputtb
                                i7 = 1;
                                %initialize i7 after last round of counting
                                while i7 < 42
                                %write into output table
                                    t = ptb(i8,i5);
                                    r = table2array(t);
                                    e = cell2mat(r);
                                    i4 = 1;
                                    while i4 < 21
                                        if e == b(i4)
                                            outputtb(i6,i7) = ptb(i8,i5);
                                            i4 = 23;
                                            %write into output table
                                        end
                                        i4 = i4 + 1;
                                    end
                                    i5 = i5 + 1;
                                    i7 = i7 + 1;
                                end
                            else
                                i10 = i10 + 1;
                                dsg(i10) = i1;
                                %pGlyco doesn't agree with trembl
                            end
                            i8 = 10000000;
                            %leave the loop of counting within trembl
                        end
                        i8 = i8 + 1;
                    end
                i4 = 1;
                end
            else
                entry1 = cell2mat(table2array(tb(i1,2)));
                %entry of peptide whose J is to the end of sequence
                prosite1 = tb(i1,3);
                prosite3 = table2array(prosite1);
                %glycosylation position
                i8 = 1;
                while i8 < (height(ptb) + 1)
                    entry2 = cell2mat(table2array(ptb(i8,5)));
                    if strcmpi(entry1,entry2)
                        gsite1 = ptb(i8,(prosite3 + 6));
                        gsite2 = table2array(gsite1);
                        gsite3 = cell2mat(gsite2);
                        if gsite3 == 'N'
                            i5 = prosite3 - 14;
                            %position of the first aa from ptb to 
                            %be written into outputtb
                            i7 = 1;
                            %initialize i7 after last round of counting
                            while i7 < 42 && i5 > 6
                            %write into output table
                                t = ptb(i8,i5);
                                r = table2array(t);
                                e = cell2mat(r);
                                i4 = 1;
                                while i4 < 21
                                    if e == b(i4)
                                        outputtb(i6,i7) = ptb(i8,i5);
                                        i4 = 23;
                                        %write into output table
                                    end
                                    i4 = i4 + 1;
                                end
                                i5 = i5 + 1;
                                i7 = i7 + 1;
                            end
                        else
                            i10 = i10 + 1;
                            dsg(i10) = i1;
                            %pGlyco doesn't agree with trembl
                        end
                        i8 = 10000000;
                        %leave the loop of counting within trembl
                    end
                    i8 = i8 + 1;
                end
                i4 = 1;
            end
            %leave write into output table loop
        end
    end
    i1 = i1 + 1;
    i6 = i6 + 1;
end