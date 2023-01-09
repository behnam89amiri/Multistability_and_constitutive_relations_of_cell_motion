function [ M ] = rectanglespec2( A, B, C, D, s)
% RECTANGLESPEC3 Creates a logical matrix of size s with ones at the
% rectangle(parallelogram) speciefied by the four points A, B, C, D. A= [a_x a_y]...


i = 1:s(2);
Mx = repmat(i,s(1),1);

i = 1:s(1);
My = repmat(i',1,s(2));

%%
%M=zeros(s);
if A(1)==B(1) % case for vertical line AB
    
%    m3=(C(2)-A(2))/(C(1)-A(1));
%    t3=C(2)-C(1)*m3;
% 
%    m4=(D(2)-B(2))/(D(1)-B(1));
%    t4=D(2)-D(1)*m4;
   
   if A(1)<C(1) 
    M= Mx>A(1) & Mx<C(1) ;  
   else
    M= Mx<A(1) & Mx>C(1) ; 
   end
       
% elseif A(1)==C(1) % case for vertical line AC
%     m1=(B(2)-A(2))/(B(1)-A(1));
%     t1=A(2)-A(1)*m1;
% 
%      m2=(D(2)-C(2))/(D(1)-C(1));
%     t2=C(2)-C(1)*m2;
%     
%     if t1<t2 && A(1)>B(1)
%     M=My-m1*Mx>t1 & My-m2*Mx<t2 & Mx<A(1) & Mx>B(1);
%     elseif t1<t2 && A(1)<B(1)
%     M=My-m1*Mx>t1 & My-m2*Mx<t2 & Mx>A(1) & Mx<B(1);    
%     elseif t1>t2 && A(1)>B(1)    
%     M=My-m1*Mx<t1 & My-m2*Mx>t2 & Mx<A(1) & Mx>B(1);    
%     else    
%     M=My-m1*Mx<t1 & My-m2*Mx>t2 & Mx>A(1) & Mx<B(1);
%     end
else % standard case with no vertical lines
    m1=(B(2)-A(2))/(B(1)-A(1));
    t1=A(2)-A(1)*m1;

     m2=(D(2)-C(2))/(D(1)-C(1));
    t2=C(2)-C(1)*m2;


%     m3=(C(2)-A(2))/(C(1)-A(1));
%     t3=C(2)-C(1)*m3;
% 
%     m4=(D(2)-B(2))/(D(1)-B(1));
%     t4=D(2)-D(1)*m4;
    
    
    if t1<t2
    M=My-m1*Mx>t1 & My-m2*Mx<t2 ;
    else
    M=My-m1*Mx<t1 & My-m2*Mx>t2 ;   
    end
end
end