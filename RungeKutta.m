%4th Order Runge Kutta method on Two Degree Freedom System

function RungeKutta

global m1 m2 c1 c2 c3 k1 k2 k3
m1 = input('Enter the value of mass1 = ')
m2 = input('Enter the value of mass2 = ')
c1 = input('Enter the value of damper coefficient1 = ')
c2 = input('Enter the value of damper coefficient2 = ')
c3 = input('Enter the value of damper coefficient3 = ')
k1 = input('Enter the value of stiffness coefficient1 = ')
k2 = input('Enter the value of stiffness coefficient2 = ')
k3 = input('Enter the value of stiffness coefficient3 = ')

t = (0:0.25:20);         % simulation time span


y0 = [-2 5 1 -3];      % [ x1(0) x2(0) x1dot(0) x2dot(0) ]   (initial values)


f = @(t,x) [x(3);               % system of 1st-order ODEs
            x(4);
            - ( (k1+k2)*x(1) - k2*x(2) + (c1+c2)*x(3) - c2*x(4) )/m1;
            - (- k2*x(1) + (k2+k3)*x(2) - c2*x(3) + (c2+c3)*x(4) )/m2];

        
yRK4 = RK4Solver(f, t, y0);     % calling Runge-Kutta 4th-order Solver
                                % plotting the numerical solution
                                
plot(t, yRK4, 'linewidth', 1.5)
grid on

%Graph Visualization
xlabel('Time, t [sec]')
ylabel({'$x_{1},\; x_{2},\; \dot{x}_{1},\; \dot{x}_{2}$'}, 'Interpreter', 'latex')
title('Time responses of the system states (with RK4 Solver)')
legend({'$x_{1}$', '$x_{2}$', '$\dot{x}_{1}$', '$\dot{x}_{2}$'}, 'Interpreter', 'latex', 'location', 'best')

end

function y = RK4Solver(f, x, y0)

y(:, 1) = y0;                    % initial condition
h = (x(2) - x(1));                 % step size
n = length(x);                   % number of steps

    for j = 1 : n-1
        
        k1 = h*f(x(j), y(:, j)) ;
        k2 = h*f(x(j) + h/2, y(:, j) + (k1/2)) ;
        k3 = h*f(x(j) + h/2, y(:, j) + (k2/2)) ;
        k4 = h*f(x(j) + h, y(:, j) + (k3/2)) ;
        k = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        y(:, j+1) = y(:, j) + (k);
    end
end

