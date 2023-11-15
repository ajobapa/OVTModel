% Example of Supercritical Hopf bifurcation
% From a stable node to a stable limit cycle shown in the phase-plane.
% Phase-plane plots of the components [x_1,x_2] from p= -2 to p=2. global p
for p = -2:0.5:2;
  tfinal=1000;
  % First set of initial values
  pts = 8;
  theta = 0:2*pi/pts:2*pi-2*pi/pts;
  x1 = 2*[cos(theta); sin(theta)];
  % Second set of initial values
  x2 = 0.01*[cos(theta); sin(theta)];
  figure(1);clf; hold on;
  for i = 1:pts
    [t x1t] = ode23(@hopf,[0 tfinal],x1(:,i));
    plot(x1t(:,1),x1t(:,2),'b')
    [t x2t]=ode23(@hopf,[0 tfinal],x2(:,i));
    plot(x2t(:,1),x2t(:,2),'r')
  end
  title(['Hopf bifurcation for p = ', num2str(p)]);
  grid on;
  xlabel('x_1'); ylabel('x_2'); box on;
  pause(1) % pause for one second so that the figure can be seen
end
