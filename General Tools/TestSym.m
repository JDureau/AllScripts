function [] = TestSym(Mat) 

if not(mean(mean(Mat == Mat'))== 1)
    disp('NotSym')
    die
end
