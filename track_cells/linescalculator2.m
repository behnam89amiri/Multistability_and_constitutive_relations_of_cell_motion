function [A] = linescalculator2(step1,step2,angle,rho,xim,yim)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% xstep1=(step1/sin(angle)+rho/sin(angle):(step1+step2)/sin(angle):xim);
% xstep2=(rho/sin(angle):(step1+step2)/sin(angle):xim);
% xsteps=sort(union(xstep1,xstep2));

%if angle<90
N=floor(sqrt(xim^2+yim^2));

rhostep1=(rho-step2-round(N/(step1+step2))*(step1+step2):(step1+step2):N);
rhostep2=(rho-round(N/(step1+step2))*(step1+step2):(step1+step2):N);
rhosteps=sort(union(rhostep1,rhostep2));% all possible distances from origin 


x1=NaN(length(rhosteps),1);
x2=NaN(length(rhosteps),1);
y1=NaN(length(rhosteps),1);
y2=NaN(length(rhosteps),1);
t_end = NaN(length(rhosteps),1);
m=-tan(angle);% slope of the line y=m*x+t with the origin on the upper left corner. Thus y values are negative and are later made positive for correct plotting of the lines
 
for i=1:length(rhosteps) % for all distances rhostep calculate the x and y points of the lines at the borders of the image
    t=rhosteps(i)*cos(angle)+rhosteps(i)*sin(angle)*tan(angle);% intersect with y axis (should be the same as rhostep(i)/cos(angle) ) 
    
    if angle==pi/2 && rhosteps(i)>0 && rhosteps(i)<xim
        x1(i)=rhosteps(i);
        x2(i)=rhosteps(i);
        y1(i)=0;
        y2(i)=yim;
    
    elseif angle==pi/2
    elseif angle==0 && rhosteps(i)>0 && rhosteps(i)<yim
        x1(i)=0;
        x2(i)=xim;
        y1(i)=rhosteps(i);
        y2(i)=rhosteps(i);
    elseif angle==0
    elseif t<=0 && t>=-yim
            x1(i)=0;
            y1(i)=-t;
        if m*xim+t>0
            x2(i)=-t/m;
            y2(i)=0;
        elseif m* xim+t<=0 && m*xim+t >=-yim
            x2(i)=xim;
            y2(i)= -m*xim-t;
        else
            x2(i)= (-yim-t)/m;
            y2(i)= yim;
        end
    elseif t>0
        if m*xim+t<=-yim
          x1(i)=-t/m;
          y1(i)=0;
          x2(i)= (-yim-t)/m;
          y2(i)= yim;
        elseif m*xim+t>-yim &&m*xim+t<0
          x1(i)=-t/m;
          y1(i)=0;  
          x2(i)=xim;
          y2(i)= -m*xim-t;
        else
        end
    elseif t<-yim
        if m*xim+t>0
          x1(i)=-t/m;
          y1(i)=0;
          x2(i)= (-yim-t)/m;
          y2(i)= yim;
        elseif m*xim+t>-yim &&m*xim+t<0
          x1(i)=(-yim-t)/m;
          y1(i)=yim;  
          x2(i)=xim;
          y2(i)= -m*xim-t;
        else
        end
        
    end
    t_end(i) = t;
end
m(1:size(x1,1),1) = m;

A=[x1 x2 y1 y2 m t_end];


%    for i=1:xim/floor(step1+step2)
%        bla=
%        
%    end
%end



end

