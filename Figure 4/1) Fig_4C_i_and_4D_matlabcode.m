%% FASTQ Read
tic

[~,S1,Q1] = fastqread('E-PCR_1.fastq');
[~,S2,Q2] = fastqread('E-PCR_2.fastq');

S1S = [S1 S2];
S1Q = [Q1 Q2];
toc

for a= 1:length(S1S)
    L1(a)= length(S1S{a});
end

clearvars -except S1S S1Q

L = zeros(length(S1S),1);
for a= 1 : length(S1S)
    L(a) = length(S1S{a});
end
histogram(L)
title('DNA length Histogram')

if exist('S1Q')
    Qsort = cell(length(S1S),1);
    for a = 1 : length(S1S)
        if length(S1S{a}) == 100
            M = zeros(1,100);
            for b = 1 : 100
                if double(S1Q{a}(b)-33) > 10
                    M(b) = 1;
                end
            end
            if sum(M) == 100
                Qsort{a} = S1S{a}(1:100);
            end
        end
    end
    empty = ~cellfun('isempty', Qsort);
    Qz = Qsort(empty);
else
    Qz = S1S;
    clear S1S S1Q
end

% Reverse DNA Align

D = 'ATGCACTCAGACGTACAACGCAGGCAGCCAAACGAACACACGCTCTGTCACCGTTCCACATTTCTATTGAGACGAGAGTACATTACACCTCACGCCTTCG';
F_primers = 'ATGCACTCAGACGTACAACG';
R_primers = 'CATTACACCTCACGCCTTCG';

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
S = cell(1,length(Qz));
parfor a = 1 : length(Qz)
    if sum(Qz{a}(end-length(rev)+1:end) == rev) >= length(rev)-2
        Reverse = Qz{a};
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
        S{a} = Qz{a};
    end
end

F = cell(length(S),1);
for a = 1 : length(S)
    W = S{a};
    if length(W) == 100
        if sum(W(1:20)==F_primers) >= 15
            if sum(W(81:100)==R_primers)>= 15
                F{a} = W;
            end
        end
    end
end
empty = ~cellfun('isempty', F);
J = F(empty);

DNAmatrix = cell2mat(J);

%최빈값 구하는 것
A = mode(cell2mat(J));

%% Codon Bar
W = DNAmatrix;
for b = 1 : 100 
    U(b,1) = length(strfind(W(:,b)','A')); % A가 포함되어 있는 개수
    U(b,2) = length(strfind(W(:,b)','T')); % T가 포함되어 있는 개수
    U(b,3) = length(strfind(W(:,b)','C')); % C가 포함되어 있는 개수
    U(b,4) = length(strfind(W(:,b)','G')); % G가 포함되어 있는 개수
end
%%
for c = 1 : length(U)
    K(c,1) = U(c,1)/sum(U(c,:))*100; % codon 비율 백분율로 표현
    K(c,2) = U(c,2)/sum(U(c,:))*100;
    K(c,3) = U(c,3)/sum(U(c,:))*100;
    K(c,4) = U(c,4)/sum(U(c,:))*100;
end
figure
bar(K,'stacked')
title("title name")
legend('A','T','C','G')
