function Qstd = stdQscore(QF,DNAlength)
Qbar1 = zeros(length(QF),100);
for a = 1 : length(QF)
    W = QF{a};
    for b = 1 : min(length(W), 100)
        Qbar1(a,b) = double(W(b)-33);
    end
end

QF1 = Qbar1(:,1:end);

Qstd1 = zeros(1,size(QF1,2));
for a = 1 : size(QF1,2)
    Qstd1(a) = std(QF1(:,a));
end

Qstd = [Qstd1]';