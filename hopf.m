
%function xdot=hopf(t,x)
%  global p
%  mu = 1; a = 1; xdot = zeros(2,1);
%  xdot(1)=-mu.*x(2)+x(1).*(p-a.*(x(1)^2)-a.*(x(2)^2));
%  xdot(2)=mu.*x(1)+x(2).*(p-a.*(x(1)^2)-a.*(x(2)^2));
%end

function xdot=hopf(x,y)
  global b
  xdot(1)=2*x-10*b*x;
  xdot(2)=10*b*x-y*(100/18)