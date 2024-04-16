% ODE/convergence parameters
tmin = 0; tmax = 1; % time interval
delta = 0.001; % tolerance parameter
N = 1000; % small time increment
t = linspace(tmin,tmax,N+1); % N+1 includes 0,

h = 1/N; % step size
h2 = h/2; % convenient for rk4

% OCP parameters (PART 7.1)
r = 1;
A = 1;
B = 12;
C = 1;
x0 = 1; % first initial condition

% initialization of x,u,lambda
u = zeros(N+1,1); % initial guess, if dividing by u, make u =/ 0
x = zeros(N+1,1); % stores solutions
lambda = zeros(N+1,1);
x(1) = x0;

% transversality condition
lambda(N+1) = C;

% initialize convergence criteria
test = -1;

% storage for iterates of u
u_iterates = [];

% - - - - - - - START SOLVING - - - - - - - for x0 = 1
while (test < 0)

% storage for previous approximations
oldu = u;
oldx = x;
oldlambda = lambda;

% solve for x (forward) using RK4 algorithm
xPrime = @(t,x,u) r.*x + A.*x.*u - B.*u.^2.*exp(-x);

for i = 1:N
u2 = 0.5*(u(i)+u(i+1)); % midpoint between u-vector values

k1 = xPrime(t(i),x(i),u(i));
k2 = xPrime(t(i)+h2,x(i)+h2*k1,u2);
k3 = xPrime(t(i)+h2,x(i)+h2*k2,u2);
k4 = xPrime(t(i)+h,x(i)+h*k3,u(i+1));

x(i+1) = x(i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

% solve for lambda (backward)
lambdaPrime = @(t,x,u,lambda) -lambda*r - lambda.*A.*u - lambda.*B.*u.^2.*exp(-x);

for i = 1:N % j index is backward movement
j = N + 2 - i;
u2 = 0.5*(u(j)+u(j-1)); % evaluating midpoints
x2 = 0.5*(x(j)+x(j-1));

k1 = lambdaPrime(t(j),x(j),u(j),lambda(j)); % subtraction for backward movement
k2 = lambdaPrime(t(j) - h2,x2, u2, lambda(j)-h2*k1);
k3 = lambdaPrime(t(j) - h2,x2, u2, lambda(j)-h2*k2);
k4 = lambdaPrime(t(j) - h, x(j-1), u(j-1),lambda(j)-h*k3);

lambda(j-1) = lambda(j) - h/6*(k1 + 2*k2 + 2*k3 + k4);
end

% update u
u1 = (lambda.*A.*x)./(2.*(1+lambda.*B.*exp(-x))); % u star
u = 0.5*(u1+oldu);

% store iterates of u
u_iterates = [u_iterates,u];

temp1 = delta*sum(abs(u)) - sum(abs(oldu - u));
temp2 = delta*sum(abs(x)) - sum(abs(oldx - x));
temp3 = delta*sum(abs(lambda)) - sum(abs(oldlambda - lambda));

test = min(temp1, min(temp2, temp3));
end

% complete solution
Y1 = [t',x,lambda,u];
