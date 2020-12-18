function Xinyi()
clear all
%% Initial condition
Nx1 = 10;
Nx2 = 20;
Lx = 15;
alpha = 1;

dx1 = Lx/(Nx1-1);
dx2 = Lx/(Nx2-1);

x1 = 0:dx1:Lx;
x2 = 0:dx2:Lx;
dt = 0.001;

% Set initial condition for T
T1_steady = zeros(1,Nx1);
T2_steady  = zeros(1,Nx2);

T1_Euler = zeros(1,Nx1);
T1_steady_check_Euler  = zeros(1,Nx1);
T2_Euler = zeros(1,Nx2);
T2_steady_check_Euler  = zeros(1,Nx2);

T1_Crank = zeros(1,Nx1);
T1_Crank_old = zeros(1,Nx1);
T1_Crank_check = zeros(1,Nx1);
T2_Crank = zeros(1,Nx2);
T2_Crank_old = zeros(1,Nx2);
T2_Crank_check = zeros(1,Nx2);

% Set exact solution
for i=1:Nx1
    T1_steady(i) = x1(i)*x1(i)*exp(-x1(i));
end

for i=1:Nx2
    T2_steady(i) = x2(i)*x2(i)*exp(-x2(i));
end

% Set boundary condition
for i=1:Nx1
T1_Euler(i,1) = 0; % T(x,0) = 0
T1_Euler(1,i) = 0; % T(0,t) = 0
end

for i=1:Nx2
T2_Euler(i,1) = 0; % T(x,0) = 0
T2_Euler(1,i) = 0; % T(0,t) = 0
end

for i=1:Nx1
T1_Crank(i,1) = 0; % T(x,0) = 0
T1_Crank(1,i) = 0; % T(0,t) = 0
T1_Crank_old(i,1) = 0; % T(x,0) = 0
T1_Crank_old(1,i) = 0; % T(0,t) = 0
end

for i=1:Nx2
T2_Crank(i,1) = 0; % T(x,0) = 0
T2_Crank(1,i) = 0; % T(0,t) = 0
T2_Crank_old(i,1) = 0; % T(x,0) = 0
T2_Crank_old(1,i) = 0; % T(0,t) = 0
end


for iteration = 1:1000
  iteration;
  
% Explicit Euler solution, Nx is 10
for i=2:Nx1
T1_Euler(i) = T1_Euler(i) + dt*(alpha*(T1_Euler(i+1)-2*T1_Euler(i)+T1_Euler(i-1))/(dx1*dx1)-(x1(i)*x1(i)-4*x1(i)+2)*exp(-x1(i)));
% Steady state checking
T1_steady_check_Euler(i) =  alpha*(T1_Euler(i+1)-2*T1_Euler(i)+T1_Euler(i-1))/(dx1*dx1)-(x1(i)*x1(i)-4*x1(i)+2)*exp(-x1(i));
end
 T1_Euler(floor(Nx1/2));

% Explicit Euler solution, Nx is 20
for i=2:Nx2
T2_Euler(i) = T2_Euler(i) + dt*(alpha*(T2_Euler(i+1)-2*T2_Euler(i)+T2_Euler(i-1))/(dx2*dx2)-(x2(i)*x2(i)-4*x2(i)+2)*exp(-x2(i)));
% Steady state checking
T2_steady_check_Euler(i) =  alpha*(T2_Euler(i+1)-2*T2_Euler(i)+T2_Euler(i-1))/(dx2*dx2)-(x2(i)*x2(i)-4*x2(i)+2)*exp(-x2(i));
end
 T2_Euler(floor(Nx2/2));

% Crank-Nicolson, Nx is 10 
for i=2:Nx1
T1_Crank(i) = T1_Crank_old(i) + dt*((alpha/2)*((T1_Crank_old(i+1)-2*T1_Crank_old(i)+T1_Crank_old(i-1)+T1_Crank(i+1)-2*T1_Crank(i)+T1_Crank(i-1))/(x1(i)*x1(i)))-(x1(i)*x1(i)-4*x1(i)+2)*exp(-x1(i)));
% Steady state checking
T1_Crank_check(i) =  (alpha/2)*((T1_Crank_old(i+1)-2*T1_Crank_old(i)+T1_Crank_old(i-1)+T1_Crank(i+1)-2*T1_Crank(i)+T1_Crank(i-1))/(x1(i)*x1(i)))-(x1(i)*x1(i)-4*x1(i)+2)*exp(-x1(i));
end
T1_Crank_old = T1_Crank;
T1_Crank(floor(Nx1/2)); 

% Crank-Nicolson, Nx is 20 
for i=2:Nx2
T2_Crank(i) = T2_Crank_old(i) + dt*((alpha/2)*((T2_Crank_old(i+1)-2*T2_Crank_old(i)+T2_Crank_old(i-1)+T2_Crank(i+1)-2*T2_Crank(i)+T2_Crank(i-1))/(x2(i)*x2(i)))-(x2(i)*x2(i)-4*x2(i)+2)*exp(-x2(i)));
% Steady state checking
T2_Crank_check(i) =  (alpha/2)*((T2_Crank_old(i+1)-2*T2_Crank_old(i)+T2_Crank_old(i-1)+T2_Crank(i+1)-2*T2_Crank(i)+T2_Crank(i-1))/(x2(i)*x2(i)))-(x2(i)*x2(i)-4*x2(i)+2)*exp(-x2(i));
end
T2_Crank_old = T2_Crank;
T2_Crank(floor(Nx2/2)); 
end

figure(1)
plot(x2,T2_steady);
xlabel('x','FontSize',12, 'FontWeight','bold','Color','k')
ylabel('T_steady','FontSize',12, 'FontWeight','bold','Color','k')

figure(2)
plot(x1,T1_Euler);hold on;
plot(x2,T2_Euler);
plot(x1,T1_steady);
plot(x2,T2_steady);
xlabel('x','FontSize',12, 'FontWeight','bold','Color','k')
ylabel('T','FontSize',12, 'FontWeight','bold','Color','k')
legend('Explict Euler, Nx = 10','Explict Euler, Nx = 20','Exact Solution, Nx = 10','Exact Solution, Nx = 20')
hold off;

figure(3)
plot(x1,T1_Crank);hold on;
plot(x2,T2_Crank);
plot(x1,T1_steady);
plot(x2,T2_steady);
xlabel('x','FontSize',12, 'FontWeight','bold','Color','k')
ylabel('T','FontSize',12, 'FontWeight','bold','Color','k')
legend('Crank-Nicolson, Nx = 10','Crank-Nicolson, Nx = 20','Exact Solution, Nx = 10','Exact Solution, Nx = 20')
