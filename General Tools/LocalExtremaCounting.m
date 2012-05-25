function n = LocalExtremaCounting(v)
% Counts the number of extremas in v

if not(size(v,1))==1
    disp('Wrong argument')
    die
end

n=0;
for i = 2:length(v)-1
    if (v(i)-v(i-1))*(v(i)-v(i+1))>0
        n=n+1;
    end
end