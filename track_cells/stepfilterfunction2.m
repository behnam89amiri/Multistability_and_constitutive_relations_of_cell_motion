function [ vec ] = stepfilterfunction2(step1,step2,N )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

vec1=zeros(N,1);
for i=1: ceil(N/(step1+step2))
vec1(round((i-1)*(step1+step2)+1):round((i-1)*(step1+step2)+step1))=1;%ones(round((i-1)*(step1+step2)+step1),1);

end
vec=vec1(1:N);


end

