%% Bayesian optimization

clear;
clear obj_PID_panda_mass;
close all;
clc;
delete(gcp('nocreate'))

FigCount = 1;

%inputs

max_degrees_of_error = 10; 

max_degrees_of_error_mass = 10;

N_explorations = 15; %number of iterations

N_exploitations = 15;

toll_qerr = max_degrees_of_error*pi/180;
toll_qerr_mass = max_degrees_of_error_mass*pi/180;
%if we go over this tollerance we say that this is unstable


%% Constants  definition  
Ts = 1e-3;         % sampling time (s). We assume measurements are available at this rate
Tsim = 3;         % simulation length (s)


professor_reference = 0; % if I want to use its reference 
%load("results_panda_PID.mat")
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
const.toll_qerr_mass = toll_qerr_mass;


bound_min_mass = 0.5;
bound_max_mass = 8 ;


%% Mass optimization

% Define objective function as a function of the optimization variables
% opt_vars_mass only, to be passed to the optimizer

% initialization of results matrices
ObjectiveTrace_mass_matrix = [];
ObjectiveMinimumTrace_mass_matrix = [];

% inicialization of mass vector
current_mass_vector = [0,0,0,0];
mass_v_index = 4;
% this vector will be updated with the estimated masses, it is used to
% input the already estimated masses in the robot_eval for successive
% estimations, after the 4 estimations it will contain the results



for ii = [6 , 4 , 2 , 1]  %indexes use to move the joint BEFORE the masses in order to observe the inertial behaviour


    % Mass
    opt_vars_mass = [];
    
    opt_vars_mass = [opt_vars_mass optimizableVariable('mass', [bound_min_mass, bound_max_mass],'Type','real')]; % mass
    % we estimate only one mass at the time

    % the OBJ function takes the variable (only one), the struct const ,
    % the index of the loop which is equivalent to the joint to move , the
    % current mass vector to store the results and utilize already
    % estimated masses
    func_mass =  @(x_vars)(obj_PID_panda_mass(x_vars,const,ii,current_mass_vector));
    initial_X = []; % array2table(rho_init, 'VariableNames', opt_var_names);
    
    clear objBO; % just to reset the function inner counter
    
    results_exploration_mass = bayesopt(func_mass, opt_vars_mass,...
        'AcquisitionFunctionName', 'lower-confidence-bound', ...
        'IsObjectiveDeterministic', true,...
        'MaxObjectiveEvaluations', N_explorations,...
        'NumCoupledConstraints',0, ...
        'NumSeedPoint',10,...
        'GPActiveSetSize', 300,...
        'PlotFcn',{@plotMinObjective,@plotObjectiveEvaluationTime});
    
    ObjectiveTrace_exploration_mass = results_exploration_mass.ObjectiveTrace;
    ObjectiveMinimumTrace_exploration_mass = results_exploration_mass.ObjectiveMinimumTrace;
    % appending to a vector in order to coller both exploration and
    % exploitations togehter


    bestPoints_mass = results_exploration_mass.XAtMinObjective; % extract best intermediate mass (table)
    intermediate_mass = bestPoints_mass.mass; % table -->scalar
    dm = intermediate_mass*0.2; % calculates a dm of 20% of the current mass
    

    % Redefining bounds for explotation with the bounds at +- dm
    opt_vars_mass = [];
    
    opt_vars_mass = [opt_vars_mass optimizableVariable('mass', [intermediate_mass-dm, intermediate_mass+dm],'Type','real')]; % mass

    
    results_mass = bayesopt(func_mass, opt_vars_mass,...
        'Verbose', 0,...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'IsObjectiveDeterministic', true,...
        'MaxObjectiveEvaluations', N_exploitations,...
        'MaxTime', inf,...
        'NumCoupledConstraints', 0, ...
        'NumSeedPoint', 10,...
        'GPActiveSetSize', 300,...
        'InitialX', bestPoints_mass,...
        'PlotFcn',{@plotMinObjective, @plotObjectiveEvaluationTime});
    
    ObjectiveTrace_mass_ii = [ObjectiveTrace_exploration_mass; results_mass.ObjectiveTrace];
    ObjectiveMinimumTrace_mass_ii = [ObjectiveMinimumTrace_exploration_mass; results_mass.ObjectiveMinimumTrace]
    % appends exploitations to exploration to form a single vector for
    % plotting

    best_vars_mass_ii = results_mass.bestPoint; % extract best point  (table)
    mass_value = best_vars_mass_ii.mass; % table --> scalar

    current_mass_vector(mass_v_index) = mass_value; % update the current mass vector
    fprintf("Best mass %d = %.4f\n", mass_v_index, mass_value); % print to screen result of optimization
    mass_v_index = mass_v_index - 1; % update the index

    %filling of results matrices (global = outside the loop)
    ObjectiveTrace_mass_matrix = [ObjectiveTrace_mass_matrix ; ObjectiveTrace_mass_ii];
    ObjectiveMinimumTrace_mass_matrix = [ObjectiveMinimumTrace_mass_matrix ; ObjectiveMinimumTrace_mass_ii];
    

end

%% end of mass optimization ( the rest of the code was not used)

%% Plot of the values of the objective function over the iterations for masses

idx_min = N_explorations + results_mass.IndexOfMinimumTrace(end);
N = N_explorations + N_exploitations;
iteration = 1:N;


figure(1)
plot(iteration, ObjectiveTrace_mass, 'b-') % Valori effettivi osservati
hold on;
plot(iteration, ObjectiveTrace_mass, 'bo', 'MarkerFaceColor', 'b') % Valori effettivi osservati
hold on;
plot(ObjectiveMinimumTrace_mass, 'r', 'LineWidth', 2) % Minimo osservato fino ad ogni iterazione
plot(iteration(idx_min), ObjectiveMinimumTrace_mass(idx_min), ...
    'MarkerEdgeColor','black','MarkerFaceColor','green', ...
    'Marker','square','MarkerSize',10); % Punto migliore in assoluto

xlabel('Iteration index $i$ (-)', 'Interpreter', 'latex');
ylabel('Performance cost  $\tilde J$ (-)', 'Interpreter', 'latex');
legend('Objective function interpolation', 'Objective function evaluation' ,'Current best point', 'Overall best point');
grid on;
title("Objective Function Value for masses estimation Over Bayesian Optimization Iterations")






%% LOOP FOR OPTIMIZATION OF 7 JOINTS

best_vars = table();

for ii = 6:-1:1
    %% Definition of the objective function

    % Define objective function as a function of the optimization variables
    % opt_vars only, to be passed to the optimizer

    func =  @(x_vars)(obj_PID_panda2(x_vars,Ts, time, n_DoFs, Robot, friction(ii), toll_qerr, q_r, dq_r, ddq_r, ii));
    initial_X = []; % array2table(rho_init, 'VariableNames', opt_var_names);

    clear objBO; % just to reset the function inner counter

    %% Perform Byayesian Optimization-EXPLORATION


    % parpool('local');
    %'UseParallel', true, ...
    fprintf("\nEXPLORATION JOINT: %.0f \n", ii);

    results_exploration{ii} = bayesopt(func, opt_vars(3*ii-2 : 3*ii),...
        'AcquisitionFunctionName', 'lower-confidence-bound', ...
        'IsObjectiveDeterministic', true,...
        'MaxObjectiveEvaluations', N_explorations,...
        'NumCoupledConstraints',0, ...
        'NumSeedPoint',10,...
        'GPActiveSetSize', 300,...
        'PlotFcn',{@plotMinObjective,@plotObjectiveEvaluationTime});


    %% Perform Byayesian Optimization-EXPLOITATION

    bestPoints = results_exploration{ii}.XAtMinObjective;
    fprintf("EXPLOITATION JOINT: %.0f \n", ii);

    results{ii} = bayesopt(func, opt_vars(3*ii-2 : 3*ii),...
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

    best_vars_i = results{ii}.bestPoint;
    best_vars = [best_vars, best_vars_i]; 

    
    
end

%% Evaluate again the objective function for the best set of parameters

% Evaluate performance of the optimal design, so we simply put inside the
% function the top parameters that we have obtained 
close all



for ii = 7:-1:1
    x_vars = best_vars{1, 3*ii-2 : 3*ii}; 
    x_vars_table = array2table(x_vars, 'VariableNames', best_vars.Properties.VariableNames(1:3));

    fprintf("\nThis is what we have seen as minimum for the joint %0.f :\n", 8 - ii)
    obj_PID_panda2(x_vars_table, Ts, time, n_DoFs, Robot, friction(ii), toll_qerr, q_r, dq_r, ddq_r, ii);
end

%% Monte-carlo simulations

% const.plot_fig = true;
% const.Tsim = 10;
% N_mc = 1;
% OBJ_VEC = zeros(N_mc,1);
% for i=1:N_mc
%     OBJ_VEC(i) = obj_PID_panda(best_vars,const);
%     disp(OBJ_VEC(i));
% end

%% Print of the PID gains

% const.Tsim = 10;
% idx_min = results.IndexOfMinimumTrace(end);
% results.XTrace(idx_min,:)
% clear objBO;
% obj_val=obj_PID_panda(results.XTrace(idx_min,:), const);

%% Plot of the values of the objective function over the iterations


for ii = 7:-1:1

    idx_min_exploration_i = results_exploration{ii}.IndexOfMinimumTrace(end);
    idx_min_exploitation_i = results{ii}.IndexOfMinimumTrace(end);
    N_exploration_i = length(results_exploration{ii}.ObjectiveTrace);
    N_exploitation_i =  length(results{ii}.ObjectiveTrace);
    iteration_i = 1:N_exploration_i + N_exploitation_i;

    if(results{ii}.ObjectiveMinimumTrace(idx_min_exploitation_i) < results_exploration{ii}.ObjectiveMinimumTrace(idx_min_exploitation_i))
        idx_min = idx_min_exploitation_i;
        exploration_flag = 0;
    else
        idx_min = idx_min_exploration_i;
        exploration_flag = 1;
    end

    figure(FigCount)
    FigCount = FigCount +1;
    plot(1:N_exploration_i, results_exploration{ii}.ObjectiveTrace, '*') %Values of J in the exploration 
    hold on
    plot((N_exploration_i+1):(N_exploration_i+N_exploitation_i), results{ii}.ObjectiveTrace, 'k*') %Values of J in the expliotation
    hold on;
    %plot(results.ObjectiveMinimumTrace, 'r', 'LineWidth', 2) % Minimo osservato fino ad ogni iterazione

    if (exploration_flag == 0)
        plot(iteration_i(idx_min + N_exploration_i), results{ii}.ObjectiveMinimumTrace(idx_min), ...
            'MarkerEdgeColor','black','MarkerFaceColor','green', ...
            'Marker','square','MarkerSize',10); % Punto migliore in assoluto
    else
        plot(iteration_i(idx_min), results_exploration{ii}.ObjectiveMinimumTrace(idx_min), ...
            'MarkerEdgeColor','black','MarkerFaceColor','green', ...
            'Marker','square','MarkerSize',10); % Punto migliore in assoluto
    end


    xlabel('Iteration index $i$ (-)', 'Interpreter', 'latex');
    ylabel('Performance cost  $\tilde J$ (-)', 'Interpreter', 'latex');
    legend('Values during exploration', 'Values during exploitation', 'Overall best point');
    grid on;
    title(sprintf("Objective Function Value Over Bayesian Optimization Iterations - Joint %d", ii));
end


%% Last computation for all the data with the best values for Ki, Kd, Kp

clear t
t = time;

Kp = diag([results.bestPoint.Kp1, results.bestPoint.Kp2, results.bestPoint.Kp3, results.bestPoint.Kp4, results.bestPoint.Kp5, results.bestPoint.Kp6, results.bestPoint.Kp7]);
Ki = diag([results.bestPoint.Ki1, results.bestPoint.Ki2, results.bestPoint.Ki3, results.bestPoint.Ki4, results.bestPoint.Ki5, results.bestPoint.Ki6, results.bestPoint.Ki7]);
Kd = diag([results.bestPoint.Kd1, results.bestPoint.Kd2, results.bestPoint.Kd3, results.bestPoint.Kd4, results.bestPoint.Kd5, results.bestPoint.Kd6, results.bestPoint.Kd7]);

% Kp = diag([results{1}.bestPoint.Kp1, results{2}.bestPoint.Kp2, results{3}.bestPoint.Kp3, results{4}.bestPoint.Kp4, results{5}.bestPoint.Kp5, results{6}.bestPoint.Kp6, results{7}.bestPoint.Kp7]);
% Ki = diag([results{1}.bestPoint.Ki1, results{2}.bestPoint.Ki2, results{3}.bestPoint.Ki3, results{4}.bestPoint.Ki4, results{5}.bestPoint.Ki5, results{6}.bestPoint.Ki6, results{7}.bestPoint.Ki7]);
% Kd = diag([results{1}.bestPoint.Kd1, results{2}.bestPoint.Kd2, results{3}.bestPoint.Kd3, results{4}.bestPoint.Kd4, results{5}.bestPoint.Kd5, results{6}.bestPoint.Kd6, results{7}.bestPoint.Kd7]);


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

figure(8)
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

figure(FigCount)
FigCount = FigCount +1;

plot(t, (q_r - q_msr) * 180/pi); 
xlabel('Time [s]');
ylabel('Position Error [deg]');
legend('Joint 1','Joint 2','Joint 3','Joint 4','Joint 5','Joint 6','Joint 7');
grid on;
title("Position Tracking Error for Each Joint")

%% Difference between Reference and Measured position

figure(FigCount)
FigCount = FigCount +1;

plot(t, q_r); 
hold on
plot(t, q_msr); 
xlabel('Time [s]');
ylabel('Joint Position [rad]');
legend({'q1','q2','q3','q4','q5','q6','q7', ...
        'qmsr1','qmsr2','qmsr3','qmsr4','qmsr5','qmsr6','qmsr7'});
grid on;
title("Joint Positions: Reference vs Measured")

%% Difference between Reference and Measured speed

figure(FigCount)
FigCount = FigCount +1;

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

figure(FigCount)
FigCount = FigCount +1;

plot(t, tau_PID); % Coppie di controllo calcolate dal PID
xlabel('Time [s]');
ylabel('Control Torque [Nm]');
legend('τ_{ctrl,1}','τ_{ctrl,2}','τ_{ctrl,3}','τ_{ctrl,4}','τ_{ctrl,5}','τ_{ctrl,6}','τ_{ctrl,7}');
grid on;
title("Control Torque from PID Controller")


%%

OBJ_VEC = [OBJ_VEC; obj_val];

save results_panda_PID results
