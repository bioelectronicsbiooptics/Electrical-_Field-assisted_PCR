function Seq = evertComplementary(Primer,Fastq)
Reverse = Primer;
A_F = strfind(Reverse,'A');
T_F = strfind(Reverse,'T');
G_F = strfind(Reverse,'G');
C_F = strfind(Reverse,'C');
S_comp = Reverse;
S_comp(A_F) = 'T';
S_comp(T_F) = 'A';
S_comp(G_F) = 'C';
S_comp(C_F) = 'G';
Complementary = flip(S_comp);

Seq = cell(1,length(Fastq));
for a = 1 : length(Fastq)
    if sum(Fastq{a}(1:length(Complementary)) == Complementary) >= length(Complementary)-2 ...
            || sum(Fastq{a}(2:1+length(Complementary)) == Complementary) >= length(Complementary)-2
            Reverse = Fastq{a};
            A_F = strfind(Reverse,'A');
            T_F = strfind(Reverse,'T');
            G_F = strfind(Reverse,'G');
            C_F = strfind(Reverse,'C');
            S_comp = Reverse;
            S_comp(A_F) = 'T';
            S_comp(T_F) = 'A';
            S_comp(G_F) = 'C';
            S_comp(C_F) = 'G';
            Seq{a} = flip(S_comp);
        else
            Seq{a} = Fastq{a};
        end
    end
end