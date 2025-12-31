function meanQ = QscoreMean(Qdata,DNAlength)
intQ = zeros(length(Qdata),100);
for a = 1 : length(Qdata)
    charData = Qdata{a};
    for b = 1 : length(charData)
        intQ(a,b) = double(charData(b)-33);
    end
end
meanQ = zeros(DNAlength,1);
for a = 1 : DNAlength
    meanQ(a) = mean(intQ(:,a));
end