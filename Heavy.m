function y = Heavy(x)
% 0 si x<0
% 1 si x>=0

if x==0
   y=1;
else
	y=(1+abs(x)/x)/2;
end;