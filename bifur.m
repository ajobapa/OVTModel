%BIFURCATION
%BIFUR(FCN,XRANGE,PRANGE) draws the bifurcation diagram for
%function FCN over specified x and p (parameter) ranges
%XRANGE is a row vector [XMIN XMAX]
%PRANGE is a row vector [PMIN PMAX]


function y = saddlefun(x,p)
  %p=1;
  %x=2;
  y=p+x.^2;
end
 
%function bifur(fcn,xrange,prange)
function bifur(fcn,xrange,prange)
  nn=100; %num of points plotted
  p1=[prange(1):(prange(2)-prange(1))/nn:prange(2)];
  x1=[xrange(1):(xrange(2)-xrange(1))/nn:xrange(2)];
  [p,x] = meshgrid(p1,x1); %generates points inbuilt function
  fval = feval(fcn,x,p);
  figure(1);
  [c,h] = contour(p,x,fval,[0,0],'r'); %plot 0 contour line
  %xlabel('p') ylabel('x')
  x=x(:);
  p=p(:);
  ind = find(abs(fval)>0.1*mean(abs(fval(:))));
  x=x(ind);
  p=p(ind);
  figure(1);
  hold on;
  plot(p,x,'go') %draw initial points
  for iter = 1:100
    x = x+0.1*feval(fcn,x,p); %solve ODE by Euler method w step size
    equals 0.1
  end
  hold on;
  plot (p(:),x(:),'bo')
end

bifur(@saddlefun, [-5 5], [-5 5]);
