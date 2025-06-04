
T=readtable('german_credit.csv');

y0 = table2array(T(:,1));

T(:,1)=[];

age = T.AgeInYears;

N=length(age);

isnotSens=zeros(N,1);

for i=1:N
    if age(i)>25
        isnotSens(i)=1;
    end
    % if T.PersonStatusSex(i)==2 || T.PersonStatusSex(i)==3
    %     isnotSens(i)=1;
    % end
end

isSens = logical(1-isnotSens);

T.AgeInYears = [];

[N,D] = size(T);

X0 = [ones(N,1),table2array(T(:,1:end-1))];

