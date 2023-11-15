%Aryav OV Treatment Model

%Algorithm 1

%Number of points in x
N = 48+1;
dx = 1/(N-1);
xx = (0:(N-1))'*dx;
dt = 0.5*(dx)^2
Tmax = 30*24/100;
max_iter = ceil(Tmax/dt)+1;
dt = Tmax/(max_iter-1);
t = (0:(max_iter-1))*dt
R = zeros(1,max_iter)
U = zeros(N,1); X = zeros(N,1); Y = zeros(N,1); V = zeros(N,1);
U1 = zeros(N,1); X1 = zeros(N,1); Y1 = zeros(N,1); V1 = zeros(N,1);
U2 = zeros(N,1); X2 = zeros(N,1); Y2 = zeros(N,1); V2 = zeros(N,1);

function x = tridiagSolve(l,d,u,b)
  n = length(d);
  for i = 2:n
    ratio = l(i-1)/d(i-1);
    d(i) = d(i) - ratio*u(i-1);
    b(i) = b(i) - ratio*b(i-1);
  end
  x = zeros(n,1);
  x(n) = b(n)/d(n);
  for i = n-1:-1:1
    x(i) = (b(i)-u(i)*x(i+1))/d(i);
  end
  x = x;
end

lambda = 2.0;
plot_iter = 0;
col = 'bgrcmybgrcmy';
%for burst = [25 50 100 150 200];
%for burst = [4 9 15 25 50];
for burst = [5 6 9 18 27];
  plot_iter = plot_iter+1;
  beta = 0.07*burst;
  D = 3.6;
  delta = 100/18;
  gamma = 2.5;
  mu = 100/48;
  theta = 1;
  
  %initialisation
  %mm at time (t) = 0
  R(1) = 2;
  
  %first iteration
  iter=1;
  X(:) = 0.84;
  Y(:) = 0.10;
  a = 5
  V(:) = (a*exp(-xx.^2/4));
  
  U(1) = 0;
  for i = 1:N-1
    Fp = (lambda * X(i+1) - mu * (theta- X(i+1) - Y(i+1) ) )/theta;
    Fm = (lambda * X(i) - mu * (theta-X(i)-Y(i)))/theta;
    U(i+1) = (xx(i)^2*U(i)+R(iter)/2*dx*(xx(i+1)^2*Fp+xx(i)^2*Fm))/(xx(i+1)^2);
  end
  R(iter+1) = R(iter)+dt*(U(N));
    
  for i = 1:N
    A = (U(i)-xx(i)*U(N))/R(iter+1);
    F = (lambda*X(i)-mu*(theta-X(i)-Y(i)))/theta;
    if i==1 || i==N
      X1(i) = X(i) + dt*(lambda*X(i)-beta*V(i)*X(i)-F*X(i));
      Y1(i) = Y(i) + dt*(beta*V(i)*X(i)-delta*Y(i)-F*Y(i));
    else
      X1(i) = X(i) - dt*A*(X(i)-X(i-1))/dx+dt*(lambda*X(i)-beta*V(i)*X(i)-F*X(i));
      Y1(i) = Y(i) - dt*A*(Y(i)-Y(i-1))/dx+dt*(beta*V(i)*X(i)-delta*Y(i)-F*Y(i));
    end
   
  end
  
  U1(1) = 0
  for i = 1:N-1
    Fp = (lambda*X1(i+1)-mu*(theta-X1(i+1)-Y1(i+1)));
    Fm = (lambda*X1(i)-mu*(theta-X1(i)-Y1(i)));
    U1(i+1) = (xx(i)^2*U1(i)+R(iter+1)/2*dx*(xx(i+1)^2*Fp+xx(i)^2*Fm))/(xx(i+1)^2);
  end
  
  A1 = -(xx.*U1(N)/R(iter+1)+2*D./(R(iter+1)^2)./xx);
  A1(1) = 0;
  A2 = -D./R(iter+1)^2;
  
  S = 2*V+2*delta*dt*Y1;
  L1 = 2*dt/dx/dx*(A2)-dt/dx*A1;
  D1 = (2-4*dt/dx^2*A2+2*gamma*dt)*ones(N,1);
  UU1 = 2*dt/dx^2*A2+dt/dx*A1;
  L1(N) = 4*dt/dx/dx*(A2)-dt/dx*A1(N);
  UU1(1) = 4*dt/dx^2*A2+dt/dx*A1(1);
  V1 = tridiagSolve(L1(2:N), D1(1:N), UU1(1:N-1), S);
  
  for iter = 2:max_iter-1
    
    if mod(iter,20)==0
      figure(1);
      plot(xx,X1,xx,Y1,xx,U1,xx,V1)
      legend('X','Y','U','V')
      title(['T = ' num2str((iter-1)*dt) ' R = ' num2str(R(iter))]);
      drawnow
    end
    
    R(iter+1) = R(iter)+0.5*dt*(3*U1(N)-U(N));
    for i = 1:N
      A = (U1(i)-xx(i)*U1(N))/R(iter+1);
      F = (lambda*X1(i)-mu*(theta-X1(i)-Y1(i)));
      if (i==1) || (i==N)
        X2(i) = X(i) + 2*dt*(lambda*X1(i)-beta*V1(i)*X1(i)-F*X1(i));
        Y2(i) = Y(i) + 2*dt*(beta*V1(i)*X1(i)-delta*Y1(i)-F*Y1(i));
        V2(i) = (4*V1(i)-V(i) + 2*dt*delta*Y2(i))/(3+2*gamma*dt);
      else
        X2(i) = X(i) - 2*dt*A*(X1(i+1)-X1(i-1))/2/dx+2*dt*(lambda*X1(i)-beta*V1(i)*X1(i)-F*X1(i));
        Y2(i) = Y(i) - 2*dt*A*(Y1(i+1)-Y1(i-1))/2/dx+2*dt*(beta*V1(i)*X1(i)-delta*Y1(i)-F*Y1(i));     
      end
    
    end
    
    U2(1) = 0;
    for i = 1:N-1
      Fp = (lambda*X2(i+1)-mu*(theta-X2(i+1)-Y2(i+1)))/theta;
      Fm = (lambda*X2(i)-mu*(theta-X2(i)-Y2(i)))/theta;
      U2(i+1) = (xx(i)^2*U2(i)+R(iter+1)/2*dx*(xx(i+1)^2*Fp+xx(i)^2*Fm))/(xx(i+1)^2);
    end
    
    A1 = -(xx.*U2(N)/R(iter+1)+2*D./(R(iter+1)^2)./xx);
    A1(1) = 0;
    A2 = -D./R(iter+1)^2;
    S = 4*V1-V+2*delta*dt*Y2;
    L1 = 2*dt/dx/dx*A2-dt/dx*A1;
    D1 = (3-4*dt/dx^2*A2+2*gamma*dt)*ones(N,1);
    UU1 = 2*dt/dx^2*A2 + dt/dx*A1;
    L1(N) = 4*dt/dx/dx*(A2)-dt/dx*A1(N);
    UU1(1) = 4*dt/dx^2*A2 + dt/dx*A1(1);
    V2 = tridiagSolve(L1(2:N), D1(1:N), UU1(1:N-1), S);
    X = X1;
    Y=Y1;
    U = U1;
    V = V1;
    X1 = X2;
    Y1 = Y2;
    U1 = U2;
    V1 = V2;
    if mod(iter,50) == 0
      X1 = 0.5*(X+X2);
      Y1 = 0.5*(Y+Y2);
      U1 = 0.5*(U+U2);
      V1 = 0.5*(V+V2);
    end
  end
  figure(2);
  hold on;
  plot(t*100/24,R,col(plot_iter));
  title('R')
  save(['burst_' num2str(burst) '.mat'])
end

