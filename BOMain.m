%% Bayesian optimization

%% Cleanup 

clear;
close all;
clc;

%%

%inputs

max_degrees_of_error = 15; 

N_explorations = 350; %number of iterations

N_exploitations = 140;


professor_reference = 0; % if I want to use its reference

toll_qerr = max_degrees_of_error*pi/180;
%if we go over this tollerance we say that this is unstable


%% Constants  definition  
Ts = 1e-3;         % sampling time (s). We assume measurements are available at this rate
Tsim = 3;         % simulation length (s)

 

%% Robot

n_DoFs = 7;

friction = [2 2 2 2 2 2 2];

Robot = panda_robot();

%% Motion Reference

t0 = 0; % [s]
tf = 4.; % [s]
time = t0:Ts:Tsim; % [s]
freq = 1; % [Hz]

q0 = [-0.7160   -0.5850    0.3504   -1.5666    0.2241   -2.1201   -2.8398];

amp = 0.3;

if (professor_reference == 1)
    for ii=1:n_DoFs
        q_r(:,ii) = q0(ii) + 10*pi/180*sin(time*2*pi*freq);
        dq_r(:,ii) = 2*pi*freq*10*pi/180*cos(time*2*pi*freq);
        ddq_r(:,ii) = -(2*pi*freq)^2*10*pi/180*sin(time*2*pi*freq);
    end
else
    for ii = 1:n_DoFs
        q_r(:,ii)   = q0(ii) + 10*pi/180 * (sin(time*2*pi*freq - pi/2) + 1);
        dq_r(:,ii)  = 2*pi*freq * 10*pi/180 * cos(time*2*pi*freq - pi/2);
        ddq_r(:,ii) = - (2*pi*freq)^2 * 10*pi/180 * sin(time*2*pi*freq - pi/2);
    end
end


clc
disp('going to optimization...');

%these are the references that we would like to follow and we add the into
%the struct r.
r.q_r = q_r;
r.dq_r = dq_r;
r.ddq_r = ddq_r;


%% Put all the constant in a structure 

const.Ts = Ts; % inner loop sampling time
const.Tsim = Tsim;
const.time = time;

const.n_DoFs = n_DoFs;

const.r = r;

const.Robot = Robot;
const.Robot_friction = friction;
const.q_0 = q0;
const.toll_qerr = toll_qerr;

%% Bound for controller's parameters 

opt_vars = [];

bound_min_Kp = 0.;
bound_max_Kp = 1000;

bound_min_Kd = 0;
bound_max_Kd = 500;

bound_min_Ki = 0.;
bound_max_Ki = 1000.;

%% Definition of optimizable varibales

% PID parameters
opt_vars = [];

%so now we create all the variables that we want to optimize
opt_vars = [opt_vars optimizableVariable('Kp1', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki1', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd1', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative
opt_vars = [opt_vars optimizableVariable('Kp2', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki2', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd2', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative
opt_vars = [opt_vars optimizableVariable('Kp3', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki3', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd3', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative
opt_vars = [opt_vars optimizableVariable('Kp4', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki4', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd4', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative
opt_vars = [opt_vars optimizableVariable('Kp5', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki5', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd5', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative
opt_vars = [opt_vars optimizableVariable('Kp6', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki6', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd6', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative
opt_vars = [opt_vars optimizableVariable('Kp7', [bound_min_Kp, bound_max_Kp],'Type','real')]; % PID proportional
opt_vars = [opt_vars optimizableVariable('Ki7', [bound_min_Ki, bound_max_Ki],'Type','real')]; % PID integral
opt_vars = [opt_vars optimizableVariable('Kd7', [bound_min_Kd, bound_max_Kd],'Type','real')]; % PID derivative

%% Definition of the objective function

% Define objective function as a function of the optimization variables
% opt_vars only, to be passed to the optimizer

func =  @(x_vars)(obj_PID_panda(x_vars,const));
initial_X = []; % array2table(rho_init, 'VariableNames', opt_var_names);

clear objBO; % just to reset the function inner counter

%% Perform Byayesian Optimization-EXPLORATION


% parpool('local');
%'UseParallel', true, ...

results_exploration = bayesopt(func, opt_vars,...
    'AcquisitionFunctionName', 'lower-confidence-bound', ...
    'IsObjectiveDeterministic', true,...
    'MaxObjectiveEvaluations', N_explorations,...
    'NumCoupledConstraints',0, ...
    'NumSeedPoint',10,...
    'GPActiveSetSize', 300,...
    'PlotFcn',{@plotMinObjective,@plotObjectiveEvaluationTime});

ObjectiveTrace_exploration = results_exploration.ObjectiveTrace;
ObjectiveMinimumTrace_exploration = results_exploration.ObjectiveMinimumTrace;


%% Perform Byayesian Optimization-EXPLOITATION

bestPoints = results_exploration.XAtMinObjective;

%new kp range according to the result of the exploration


bound_min_Kp_exploit = bestPoints{1, 1} - bestPoints{1, 1}*(15/100);
bound_max_Kp_exploit = bestPoints{1, 1} + bestPoints{1, 1}*(15/100);

%new ki range according to the result of the exploration
bound_min_Ki_exploit = bestPoints{1, 2} - bestPoints{1, 2}*(15/100);
bound_max_Ki_exploit = bestPoints{1, 2} + bestPoints{1, 2}*(15/100);

%new kd range according to the result of the exploration
bound_min_Kd_exploit = bestPoints{1, 3} - bestPoints{1, 3}*(15/100);
bound_max_Kd_exploit = bestPoints{1, 3} + bestPoints{1, 3}*(15/100);

%if the new min goes under zero
if(bound_min_Kp_exploit < 0)
    bound_min_Kp_exploit = bound_min_Kp;
end

if(bound_min_Ki_exploit < 0)
    bound_min_Ki_exploit = bound_min_Ki;
end

if(bound_min_Kd_exploit < 0)
    bound_min_Kd_exploit = bound_min_Kd;
end

%if the new max goes above the previous limit
if(bound_max_Kp_exploit > bound_max_Kp)
    bound_max_Kp_exploit = bound_max_Kp;
end

if(bound_max_Ki_exploit > bound_max_Ki)
    bound_max_Ki_exploit = bound_max_Ki;
end

if(bound_max_Kd_exploit > bound_max_Kd)
    bound_max_Kd_exploit = bound_max_Kd;
end

opt_vars_exploit = [
    optimizableVariable(sprintf('kp%d', ii), [bound_min_Kp_exploit, bound_max_Kp_exploit],'Type','real')
    optimizableVariable(sprintf('ki%d', ii), [bound_min_Ki_exploit, bound_max_Ki_exploit],'Type','real')
    optimizableVariable(sprintf('kd%d', ii), [bound_min_Kd_exploit, bound_max_Kd_exploit],'Type','real')
    ];

results = bayesopt(func, opt_vars,...
    'Verbose', 0,...
    'AcquisitionFunctionName', 'expected-improvement', ...
    'IsObjectiveDeterministic', true,...
    'MaxObjectiveEvaluations', N_exploitations,...
    'MaxTime', inf,...
    'NumCoupledConstraints', 0, ...
    'NumSeedPoint', 10,...
    'GPActiveSetSize', 300,...
    'InitialX', bestPoints,...
    'PlotFcn',{@plotMinObjective, @plotObjectiveEvaluationTime});

ObjectiveTrace = [ObjectiveTrace_exploration; results.ObjectiveTrace];
ObjectiveMinimumTrace = [ObjectiveMinimumTrace_exploration; results.ObjectiveMinimumTrace]

best_vars = results.bestPoint; 
%correct!! because the first iteration of the exploitation starts with the
%best of the exploration so in best_vars we have also the best of the
%exploration

%% Evaluate again the objective function for the best set of parameters

% Evaluate performance of the optimal design, so we simply put inside the
% function the top parameters that we have obtained 
close all
fprintf("\nThis is what we have seen as minimum:\n")
obj_PID_panda(best_vars,const);


%% Monte-carlo simulations

const.plot_fig = true;
const.Tsim = 10;
N_mc = 1;
OBJ_VEC = zeros(N_mc,1);
for i=1:N_mc
    OBJ_VEC(i) = obj_PID_panda(best_vars,const);
    disp(OBJ_VEC(i));
end

%% Print of the PID gains

const.Tsim = 10;
idx_min = results.IndexOfMinimumTrace(end);
results.XTrace(idx_min,:)
clear objBO;
obj_val=obj_PID_panda(results.XTrace(idx_min,:), const);

%% Plot of the values of the objective function over the iterations

idx_min = N_explorations + results.IndexOfMinimumTrace(end);
N = N_explorations + N_exploitations;
iteration = 1:N;

figure(1)

[~, idx_min] = min(ObjectiveTrace);
N = N_explorations + N_exploitations;
iteration = 1:N;

hold on;

% Split the plotting  for exploration % exploitation
plot(iteration(1:N_explorations+1), ObjectiveTrace(1:N_explorations+1), 'Color', [0.5 0 0.5], 'LineStyle', '-', ...
    'DisplayName', 'Exploration (purple)');
plot(iteration(1:N_explorations+1), ObjectiveTrace(1:N_explorations+1), 'o', 'Color', [0.5 0 0.5], ...
    'MarkerFaceColor', [0.5 0 0.5],'HandleVisibility', 'off');

plot(iteration(N_explorations+1:N), ObjectiveTrace(N_explorations+1:N), 'b-', ...
    'DisplayName', 'Exploitation (blue)');
plot(iteration(N_explorations+1:N), ObjectiveTrace(N_explorations+1:N), 'bo', 'MarkerFaceColor', 'b','HandleVisibility', 'off');

% Plot the minimum trace
plot(iteration, ObjectiveMinimumTrace, 'r', 'LineWidth', 2, ...
    'DisplayName', 'Current best point');

% Highlight global best
plot(iteration(idx_min), ObjectiveMinimumTrace(idx_min), ...
    'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Overall best point');

xlabel('Iteration index $i$ (-)', 'Interpreter', 'latex');
ylabel('Performance cost  $\tilde J$ (-)', 'Interpreter', 'latex');
legend('Location', 'best');
grid on;


title(sprintf("Objective Function Value Over Bayesian Optimization Iterations"));


%% Last computation for all the data with the best values for Ki, Kd, Kp

clear t
t = time;

Kp = diag([results.bestPoint.Kp1, results.bestPoint.Kp2, results.bestPoint.Kp3, results.bestPoint.Kp4, results.bestPoint.Kp5, results.bestPoint.Kp6, results.bestPoint.Kp7]);
Ki = diag([results.bestPoint.Ki1, results.bestPoint.Ki2, results.bestPoint.Ki3, results.bestPoint.Ki4, results.bestPoint.Ki5, results.bestPoint.Ki6, results.bestPoint.Ki7]);
Kd = diag([results.bestPoint.Kd1, results.bestPoint.Kd2, results.bestPoint.Kd3, results.bestPoint.Kd4, results.bestPoint.Kd5, results.bestPoint.Kd6, results.bestPoint.Kd7]);

B               = zeros(n_DoFs,n_DoFs,length(t));
g               = zeros(length(t),n_DoFs);
tau_l           = zeros(length(t),n_DoFs);
q_msr           = zeros(length(t),n_DoFs);
dq_msr          = zeros(length(t),n_DoFs);
ddq_msr         = zeros(length(t),n_DoFs);
tau_PID         = zeros(length(t),n_DoFs);
tau_comp        = zeros(length(t),n_DoFs);

q_msr(1,:) = q_r(1,:);

B(:,:,1)  = Robot.inertia(q_r(1,:));
g(1,:)    = Robot.gravload(q_r(1,:));

tau_l(1,:) = g(1,:)' + friction*dq_msr(1,:)';

q_err = 0;

wb = waitbar(0,'Please wait...');

ierr_m(1,:) = [0 0 0 0 0 0 0];

for jj=2:length(t)
    
    B(:,:,jj)  = Robot.inertia(q_msr(jj-1,:));
    g(jj,:)    = Robot.gravload(q_msr(jj-1,:));
    
    tau_l(jj,:)    = g(jj,:)' + friction*dq_msr(jj,:)';
    tau_PID(jj,:)  = B(:,:,jj)*(Kp * (q_r(jj,:) - q_msr(jj-1,:))' - Kd * dq_msr(jj-1,:)' + Ki * ierr_m(jj-1,:)');
    tau_comp(jj,:) = g(jj,:)' + friction*dq_msr(jj,:)';
    
    Beq = B(:,:,jj);
    
    ddq_msr(jj,:) = (Beq)\( -tau_l(jj,:)' + tau_PID(jj,:)' + tau_comp(jj,:)' );
    dq_msr(jj,:) = dq_msr(jj-1,:) + ddq_msr(jj,:)*Ts;
    q_msr(jj,:) = q_msr(jj-1,:) + dq_msr(jj,:)*Ts;
    
    ierr_m(jj,:) = ierr_m(jj-1,:) + (q_r(jj,:) - q_msr(jj,:)) * Ts;
    
    q_err = q_err + abs(q_r(jj,:) - q_msr(jj,:));
    
    waitbar(jj/length(t),wb);
    
end

J_err = sum(q_err)/length(t);

close(wb)

step = 100;


%% Robot animation

figure(2)
subplot(2, 1, 1)
title("Robot Animation During Expected-Trajectory Tracking")
for jj = 1:step:length(t)
    plot(Robot, q_r(jj,:)); 
end
subplot(2, 1, 2)
title("Robot Animation During Measured-Trajectory Tracking")
for jj = 1:step:length(t)
    plot(Robot, q_msr(jj,:)); 
end

%% Position error for each joint

figure(3)
plot(t, (q_r - q_msr) * 180/pi); 
xlabel('Time [s]');
ylabel('Position Error [deg]');
legend('Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6','Joint 7');
grid on;
title("Position Tracking Error for Each Joint")

%% Difference between Reference and Measured position

figure(4)
plot(t, q_r); 
hold on
plot(t, q_msr); 
xlabel('Time [s]');
ylabel('Joint Position [rad]');
legend({'q1','q2','q3','q4','q5','q6','q7', ...
        'qmsr1','qmsr2','qmsr3','qmsr4','qmsr5','qmsr6','qmsr7'});
title("Joint Positions: Reference vs Measured")
grid on;

%% Difference between Reference and Measured speed

figure(5)
plot(t, dq_r); % Velocità desiderata
hold on
plot(t, dq_msr); % Velocità effettiva misurata
xlabel('Time [s]');
ylabel('Joint Velocity [rad/s]');
legend({'dq1','dq2','dq3','dq4','dq5','dq6','dq7', ...
        'dqmsr1','dqmsr2','dqmsr3','dqmsr4','dqmsr5','dqmsr6','dqmsr7'});
grid on;
title("Joint Velocities: Reference vs Measured")

%% Torque from PID controller

figure(6)
plot(t, tau_PID); % Coppie di controllo calcolate dal PID
xlabel('Time [s]');
ylabel('Control Torque [Nm]');
legend('τ_{ctrl,1}','τ_{ctrl,2}','τ_{ctrl,3}','τ_{ctrl,4}','τ_{ctrl,5}','τ_{ctrl,6}','τ_{ctrl,7}');
grid on;
title("Control Torque from PID Controller")

%%

OBJ_VEC = [OBJ_VEC; obj_val];

save results_panda_PID results