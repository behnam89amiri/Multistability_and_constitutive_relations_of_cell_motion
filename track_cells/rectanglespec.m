function [ M ] = rectanglespec( A, B, C, D, s)
% RECTANGLESPEC Creates a logical matrix of size s with ones at the
% rectangle(parallelogram) speciefied by the four points A, B, C, D. A= [a_x a_y]...


i = 1:s(2);
Mx = repmat(i,s(1),1);

i = 1:s(1);
My = repmat(i',1,s(2));

%%
if A(1)==B(1) % case for vertical line AB

   
   if A(1)<C(1) 
    M= Mx>A(1) & Mx<C(1) ;  
   else
    M= Mx<A(1) & Mx>C(1) ; 
   end
   
else % standard case with no vertical lines
    m1=(B(2)-A(2))/(B(1)-A(1));
    t1=A(2)-A(1)*m1;

     m2=(D(2)-C(2))/(D(1)-C(1));
    t2=C(2)-C(1)*m2;
    
    if t1<t2
    M=My-m1*Mx>t1 & My-m2*Mx<t2 ;
    else
    M=My-m1*Mx<t1 & My-m2*Mx>t2 ;   
    end
end
end