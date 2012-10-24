function [acrate]=acrate_ov(theta);
acrate=zeros(1,size(theta,2));
for i=1:size(theta,2)
    k1=1*(diff(theta(:,i))~=0);
    acrate(i)=mean(k1);
end;