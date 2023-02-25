% x coordinates for ceil()
x1 = [-2 -2 -1 -1 0 0 1 1];
% y coordinates for ceil()
y1 = [-3 -2 -2 0 0 2 2 3];

% x coordinates for round()
x2 = [-3 -2.5 -2.5 -1.5 -1.5 -.5 -.5 .5 .5 1.5 1.5 2.5 2.5 3];
% y coordinates for round()
y2 = [-3 -3 -2 -2 -1 -1 0 0 1 1 2 2 3 3];

% x coordinates for floor()
x3 = [-3 -2 -2 -1 -1 0 0 1 1 2 2 3 3];
% y coordinates for floor()
y3 = [-1.5 -1.5 -1 -1 -.5 -.5 0 0 .5 .5 1 1 1.5];

% plot ceil() with a dashed line
plot(x1, y1, 'r--');
hold on;

% plot round() with a solid line
plot(x2, y2, 'm-');

% plot floor() with a dashed line with dots
plot(x3, y3, '-. ', 'Color', 'black');

% add title, axis labels, and legend
title('Graph of ceil(), round(), and floor() Functions');
xlabel('x');
ylabel('f(x)');
% Add a grid to the plot
grid on;
legend('ceil()', 'round()', 'floor()');

% Define the x values
x4 = 0:0.1:6;

% Define the y values for the first function
y4 = cos(x4);

% Create a new figure
figure;

% Plot the first function with a solid blue line
plot(x4, y4, 'm:', 'LineWidth', 2);

% Add a title to the plot
title('cos(x)');

% Add a label to the x-axis
xlabel('x');

% Add a label to the y-axis
ylabel('f(x)');

% Define the x values
x4 = 0:0.1:6;

% Define the y values for the second function
y5 = sin(x4);

% Create a new figure
figure;

% Plot the second function with a dashed red line
plot(x4, y5, 'r--', 'LineWidth', 2);

% Add a title to the plot
title('sin(x)');

% Add a label to the x-axis
xlabel('x');

% Add a label to the y-axis
ylabel('f(x)');

%Solution of non-linear equations
initial_points = [0.45, -0.1, 2.02];

for i = 1:length(initial_points)
    % Define the function handle for fsolve
    my_eq = @(x) (15*x - 6) / (x - 1);
    
    % Call fsolve with initial point and function handle
    sol = fsolve(my_eq, initial_points(i));
    
    % Display initial point and solution
    fprintf("Initial point: %f, Solution: %f\n", initial_points(i), sol);
end

% Define the equation as a function
eqn1 = @(k2) k2 - (13.95 - (13.90 - k2 + (10)*(k2^0.3)/(1.05))/(1 + 0.85*((3)*(k2^(-0.07)))));

% Set the options for the solver
options = optimoptions('fsolve','Display','off','MaxFunctionEvaluations',10000);

% Define the range of possible solutions for k2
k2_range = linspace(0.1, 100, 1000);

% Initialize vectors to store the solutions for k2, c1, and ca1
k2_solutions = zeros(length(k2_range), 1);
c1_solutions = zeros(length(k2_range), 1);
ca1_solutions = zeros(length(k2_range), 1);

% Iterate through the range of possible solutions for k2
for i = 1:length(k2_range)
    % Solve for k2 using fsolve
    k2_solutions(i) = fsolve(eqn1, k2_range(i), options);

    % Define the equation for c1 using the solution for k2
    eqn2 = @(c1) c1 - (10 * 3^0.3 + 3.4 - k2_solutions(i));

    % Solve for c1 using fsolve
    c1_solutions(i) = fsolve(eqn2, 1, options);

    % Define the equation for ca1 using the solutions for k2 and c1
    eqn3 = @(ca1) ca1 - (10* 3^0.3 - c1_solutions(i));

    % Solve for ca1 using fsolve
    ca1_solutions(i) = fsolve(eqn3, 1, options);
end

% Print the solution variables to the console
disp(['The solution for k2 is: ', num2str(k2_solutions(i))]);
disp(['The solution for c1 is: ', num2str(c1_solutions(i))]);
disp(['The solution for ca1 is: ', num2str(ca1_solutions(i))]);

% Plot the solutions for k2, c1, and ca1 as a function of the interest rate
interest_rates = linspace(0.01, 0.2, length(k2_range));
figure
plot(interest_rates, k2_solutions, interest_rates, c1_solutions, interest_rates, ca1_solutions)
xlabel('Interest rate')
legend('k2', 'c1', 'ca1')
