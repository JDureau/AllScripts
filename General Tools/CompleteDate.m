function Date = CompleteDate(Date)

MonthsInLetters = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

if not(isfield(Date,'MonthInLetters'))
    Date.MonthInLetters = MonthsInLetters{Date.Month};
elseif isempty(Date.MonthInLetters)
    Date.MonthInLetters = MonthsInLetters{Date.Month};
elseif not(isfield(Date,'Month'))
    Date.MonthInLetters = find(Date.MonthInLetters,MonthsInLetters);
end

