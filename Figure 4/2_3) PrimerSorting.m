function [SortingSeq,Qdata] = PrimerSorting(Fastq,Qscore,Primer_1,Primer2)
tempSeq = cell(length(Fastq),1);
for i = 1:length(Fastq)
    temp = Fastq{i};
    if length(temp) == 100 
        if sum(temp(1:20) == Primer_1) == 20
            if sum(temp(81:100) == Primer2) == 20
                tempSeq{i} = temp(1:100);
            end
        end
    end
end
empty = ~cellfun('isempty', tempSeq);
SortingSeq = tempSeq(empty);
Qdata = Qscore(empty);