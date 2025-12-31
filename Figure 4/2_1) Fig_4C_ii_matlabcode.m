%% Fastq read
DNAlength = input("DNA length?:");
%% Fastq read
[~,S1,Q1] = fastqread('E-PCR_1.fastq');
[~,S2,Q2] = fastqread('E-PCR_2.fastq');

S1S = [S1 S2]; S1Q = [Q1 Q2];

F1 = 'ATGCACTCAGACGTACAACG'; R1 = 'CATTACACCTCACGCCTTCG';

%% Evert Complementary Sequence
S_1 = evertComplementary(F1,S1S);
S_2 = evertComplementary(R1,S_1);
%% Primer Sorting
[ConS_1,ConQ_1] = PrimerSorting(S_2,S1Q,F1,R1);
%% Qscore Mean
meanQ_1 = QscoreMean(ConQ_1,DNAlength);
csvwrite('100bp Qscore Average.csv',meanQ_1)
%% Qscore std
Qstd = stdQscore(ConQ_1,DNAlength);
csvwrite('Qscore Standard deviation.csv',Qstd)
