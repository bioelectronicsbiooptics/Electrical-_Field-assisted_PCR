[~,S1,~] = fastqread('E-PCR_1.fastq');
[~,S2,~] = fastqread('E-PCR_2.fastq');
S1S = [S1 S2];

% template primer sequence
F_primers = 'ATGCACTCAGACGTACAACG';
R_primers = 'CATTACACCTCACGCCTTCG';

% template sequence
D_seq ='ATGCACTCAGACGTACAACGCAGGCAGCCAAACGAACACACGCTCTGTCACCGTTCCACATTTCTATTGAGACGAGAGTACATTACACCTCACGCCTTCG';
%% 변수
original = D_seq;
Seq_length = 100;
%%
Reverse = F_primers;
A_F = strfind(Reverse,'A');
T_F = strfind(Reverse,'T');
G_F = strfind(Reverse,'G');
C_F = strfind(Reverse,'C');
S_comp = Reverse;
S_comp(A_F) = 'T';
S_comp(T_F) = 'A';
S_comp(G_F) = 'C';
S_comp(C_F) = 'G';
rev = flip(S_comp);

% Complementary sequence
i = 0;
S = cell(1,length(S1S));
for a = 1 : length(S1S)
    if sum(S1S{a}(end-length(rev)+1:end) == rev) >= length(rev)-2
        Reverse = S1S{a};
        A_F = strfind(Reverse,'A');
        T_F = strfind(Reverse,'T');
        G_F = strfind(Reverse,'G');
        C_F = strfind(Reverse,'C');
        S_comp = Reverse;
        S_comp(A_F) = 'T';
        S_comp(T_F) = 'A';
        S_comp(G_F) = 'C';
        S_comp(C_F) = 'G';
        S{a} = flip(S_comp);
         i = i + 1;
    else
        S{a} = S1S{a};
    end
end

%% primer 위치별로 insertion/deletion 찾기
B = zeros(length(S),1);
for a = 1 : length(S)
    if length(S{a}) > Seq_length
        if sum(S{a}(Seq_length-18:Seq_length+1) == R_primers) >= 20
            B(a) = 1;
        end
    end
end
"Insertion +1: " +sum(B)

C = zeros(length(S),1);
for a = 1 : length(S)
    if length(S{a}) > Seq_length+1
        if sum(S{a}(Seq_length-17:Seq_length+2) == R_primers) >= 20
            C(a) = 1;
        end
    end
end
"Insertion +2: " +sum(C)

D = zeros(length(S),1);
for a = 1 : length(S)
    if length(S{a}) == Seq_length - 1
        if sum(S{a}(Seq_length-20:Seq_length-1) == R_primers) >= 20
            D(a) = 1;
        end
    end
end
"Deletion -1: " +sum(D)

E = zeros(length(S),1);
for a = 1 : length(S)
    if length(S{a}) == Seq_length - 2
        if sum(S{a}(Seq_length-21:Seq_length-2) == R_primers) >= 20
            E(a) = 1;
        end
    end
end
"Deletion -2: " +sum(E)

F = zeros(length(S),1);
for a = 1 : length(S)
    if length(S{a}) >= Seq_length
        if sum(S{a}(1:Seq_length) == original) == Seq_length - 1
            F(a) = 1;
        end
    end
end
"Substitution 1: "+sum(F)

G = zeros(length(S),1);
for a = 1 : length(S)
    if length(S{a}) >= Seq_length
        if sum(S{a}(1:Seq_length) == original) == Seq_length - 2
            G(a) = 1;
        end
    end
end
"Substitution 2: "+sum(G)