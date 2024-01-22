
%% Numerical Solver for a Flat Plate Laminar Boundary Layer
% Solve the reduced Navier Stokes equations for the
% boundary layer profile over a 10cm flat plate in air with an air
% speed of 1 m/s and compare the the Blasius analytical solution.
 clearvars
 close all

% Air properties (Standard Temperature and Pressure Conditions)
Uinf = 1.0;         % Freestream wind speed [m/s]
mu = 1.8e-5;        % Dynamic viscosity [kg/m-s]
rho = 1.225;        % Density [kg/m^3]
nu = mu/rho;        % Kinematic viscosity [m^2/s]

%Set up the domain
L = 0.1;            % Plate length (10 cm) [m]
H = L/10;           % Aribitrarily define domain height [m]
ix = 1000;          % Number of grid point in the x-direction
jx = 500;           % Number of grid point in the y-direction

% Sub-divide the domain into blocks, based on ix and jx, where the Navier
% Stokes equations are computed.
dx = L/(ix-1);      
dy = H/(jx-1);

% x- and y-points in the domain
x = 0:dx:L;
y = 0:dy:H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical solution (DO NOT TOUCH THIS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Intialize velocity matrices and apply the initial conditions
%[u = U_inf and v = 0 everywhere]
u_n = ones(ix,jx);      % u(y) at current time step, n
u_n = Uinf*u_n;
u_np1 = ones(ix,jx);    % u(y) at next time step, n+1
u_np1 = Uinf*u_np1;
v_n = zeros(ix,jx);     % v(y) at current time step, n
v_np1 = zeros(ix,jx);   % v(y) at next time step, n+1

%Apply the boundary Conditions
%no-slip at plate [u(x,y=0) = 0]
%v = 0 at the plate is already satisfied by the initial condition
u_n(:,1) = 0.0;
u_np1(:,1) = 0.0;

%%%%%%%%%%%%%%% This is the actual main computation loop %%%%%%%%%%%%% 
% Iterate until convergence achieved
maxresidual = 1e-06;
converged = 1;
n = 0;  % Iteration counter
while converged >= maxresidual
    n = n+1;
    for i = 1:ix-1
        for j = 2:jx-1
            % Compute alpha and beta
            alpha = nu*dx/(u_n(i,j)*dy^2);
            beta = (v_n(i,j)*dx)/(2*u_n(i,j)*dy);
        
            % Compute u in the next cell (i+1) at the next time step (n+1)
            u_np1(i+1,j) = (1/(1+2*alpha))*(u_n(i,j) + alpha*(u_n(i+1,j+1)+...
                u_n(i+1,j-1)) - beta*(u_n(i,j+1) - u_n(i,j-1))); 
        end
    end
    for i = 1:ix-1
        for j = 2:jx-1
            % Compute v in the next cell (i+1) at the next time step (n+1)
            v_np1(i+1,j) = v_n(i+1,j-1) - (dy/(2*dx))*(u_np1(i+1,j) - u_np1(i,j) +...
                u_np1(i+1,j-1) - u_np1(i,j-1));
        end
    end
    %Compute residual value at this time step
    A = abs(u_n - u_np1);   % How much did the solution change between time steps?
    converged = max(A(:));
    residual(n) = converged;
    %Update u and v matrices in time (solution here becomes initial
    %solution for next time step).
    u_n = u_np1; 
    v_n = v_np1;
end
u = u_n;
v = v_n;

% Here's a pretty plot solution of the solution over the whole domain
% (looks like a fancy CFD picture) to help you visualize things.
figure(1)
imagesc(x,y,u')
set(gca,'Ydir','normal')
colormap('jet')
c = colorbar;
c.Label.String = 'u(y) [m/s]';
xlabel('x [m]')
ylabel('y [m]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the boundary layer height, delta, along the length of the plate
% from the numerical solution
delta = zeros(1,ix);    % Vector of delta's for the numerical solution
for i = 1:ix    % Move along the plate in the x-direction
    jDelta = 0;
    for j = 1:jx    % Check u-velocity at each height for this x-location
        utemp = u(i,j);
        if utemp <= 0.99*Uinf   % Is u = 0.99*Uinf yet?
            jDelta = j;         % Save the index of the height
        end
    end
    % Use jDelta to find height where u = 0.99Uinf at this x-location
    delta(i) = y(jDelta);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% PART 1) Analytical solution - Blasius Boundary Layer Thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate Blasius relation for boundary layer thickness, delta,
% along the whole plate length, from x = 0 to L.

Rex=(rho*Uinf*x)/mu;
deltaB=(5.0*x)./sqrt(Rex);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2) Plot boundary layer height along the plate for the Blasius
% solution and the numerical solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(x,deltaB)
xlabel('x-values')
ylabel('Delta-values')
hold on 
plot(x,delta)
legend('deltaB','deltaN')
title('delta(x) vs x for Blasius and numerical calculations')

% Uncomment the line above and put the plot commands here :)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Drag over the flat plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: skip the first x-location at the leading-edge of the plate
% because the boundary layer has zero height there; it gives us a
% numerical error!

% PART 3) Blasius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
cf=0.664./sqrt(Rex);


twb=0.5*rho*(Uinf).^2*cf(1,2:end);

DragB=sum(twb.*dx);

% PART 4) Numerical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

twn=mu*((u(2:end,2)-u(2:end,1))/dy);

DragN=sum(twn*dx);

compare=DragN-DragB;

% Add commands to print the drag values:

fprintf("DragB=%4.6f\n",DragB)
fprintf("DragN=%4.6f\n",DragN)
fprintf("comparison=%4.6f\n",compare)
fprintf(['This is the difference between the numerical drag and the\n' ...
    'Blasius drag, as the x-values increase the difference between\n' ...
    'the deltas for Blasius and numerical increase slightly\n'])