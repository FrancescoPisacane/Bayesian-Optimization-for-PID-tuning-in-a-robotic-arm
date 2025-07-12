%function J = objMPC(x_var,const)
%
%   Evaluates the objective function to be minimized using bayesian
%   optimization
%   Inputs: 
%       X_VAR - variables to be optimized (in a table-like
%       data structure)
%       CONST - structure with constant values 

function J = obj_PID_panda(x_var, const)

    
    %% Count the number of function evaluations %%
    persistent idx_sim;
    if isempty(idx_sim)
        idx_sim = 1;
    end
    %% assign all constant values %%
    
    sim_var = const;
    % assign bayesian optimization variables (overwrite default constants
    % if there exist a field with the same name )
    assert(size(x_var,1) == 1)

    x_var_struct = table2struct(x_var(1,:));
    fname = fieldnames(x_var_struct);
    for idx_fn = 1:length(fname)
       sim_var.(fname{idx_fn}) = x_var_struct.(fname{idx_fn});
    end
    
    %we are creating a struct with all the constants and optimizable
    %variabels 

    %% Simulation settings - inner loop and outer loop configuration %%
    %% 
    % Setup PID object (inner loop controller)
    Kp = diag([sim_var.Kp1, sim_var.Kp2, sim_var.Kp3, sim_var.Kp4, sim_var.Kp5, sim_var.Kp6, sim_var.Kp7]);
    Ki = diag([sim_var.Ki1, sim_var.Ki2, sim_var.Ki3, sim_var.Ki4, sim_var.Ki5, sim_var.Ki6, sim_var.Ki7]);
    Kd = diag([sim_var.Kd1, sim_var.Kd2, sim_var.Kd3, sim_var.Kd4, sim_var.Kd5, sim_var.Kd6, sim_var.Kd7,]);
    
    Robot_eval = panda_robot();
    
    %%
    
    toll_qerr = const.toll_qerr;
    Ts = const.Ts;
    t = const.time;
    n_DoFs = const.n_DoFs;
    Robot = const.Robot;
    q_r = const.r.q_r;
    dq_r = const.r.dq_r;
    ddq_r = const.r.ddq_r;
    f = const.Robot_friction;
    
    %we extract all the constants (again) from the struct given in input
        
    %%
    % Run the simulation
    
    %we have a value for each parameter for each time instant
    g_eval = zeros(length(t),n_DoFs); 
    %this is a 2D matrix in which we will have the torque to stand still

    B_eval = zeros(n_DoFs,n_DoFs,length(t)); 
    %this inertia matrix is a 3D matrix with a squared matrix n_DoFs x n_DoFs
    %but we update it at every time instant and so we have a 3D matrix at
    %the end
    
    B            = zeros(n_DoFs,n_DoFs,length(t));
    g            = zeros(length(t),n_DoFs);
    tau_l        = zeros(length(t),n_DoFs); %feedback linearization
    q_msr        = zeros(length(t),n_DoFs); %measured position
    qerr         = zeros(length(t),n_DoFs); %error on position
    dqerr        = zeros(length(t),n_DoFs); %measured speed
    ddqerr       = zeros(length(t),n_DoFs); %error on speed
    dq_msr       = zeros(length(t),n_DoFs); %measured acceleration
    ddq_msr      = zeros(length(t),n_DoFs); %error on acceleration
    tau_PID      = zeros(length(t),n_DoFs); %torque from PID controller
    tau_comp     = zeros(length(t),n_DoFs); %torque for compensation the inertia of the robot
    
    %now we inizialize the measured position to the inizial position
    %defined in the B0Main
    q_msr(1,:) = q_r(1,:);
    
    %now we inizialize the rotational inertia(that doesn't depend on the
    %gravity) and the gravitational contribution to the initial values
    B(:,:,1)  = Robot.inertia(q_r(1,:));
    g(1,:)    = Robot.gravload(q_r(1,:));
    
    B_eval(:,:,1) = Robot_eval.inertia(q_r(1,:));
    g_eval(1,:)   = Robot_eval.gravload(q_r(1,:));
    
    %this creates a bar that tells us the process for the simulation
    wb = waitbar(0,'Please wait...');
    
    %integral error, this is the one that we have to initialize because it
    %has a "past history" that we have to consider
    ierr_m(1,:) = [0 0 0 0 0 0 0];

    jj = 1; %while loop counter
    exit_flag = false;

    penalty = [0 0 0 0 0 0 0]; %penalty for instabilities for each joint
    
    p = 0.01; %final scaling factor for log-penalty

    J = 0; %initialize the cost function

    while ( jj<length(t) && exit_flag==false)
        
        jj=jj+1;
        
       
        B(:,:,jj)  = Robot.inertia(q_msr(jj-1,:));
        g(jj,:)    = Robot.gravload(q_msr(jj-1,:));
        
        B_eval(:,:,jj) = Robot_eval.inertia(q_msr(jj-1,:));
        g_eval(jj,:)   = Robot_eval.gravload(q_msr(jj-1,:));
        
        %the torque for the feedback linearization (for standing still)
        tau_l(jj,:)    = f*dq_msr(jj-1,:)' + g(jj,:)' ;
        
        %the torque from PID controller
        tau_PID(jj,:)  = B_eval(:,:,jj)*(Kp * (q_r(jj,:) - q_msr(jj-1,:))' - Kd * dq_msr(jj-1,:)' + Ki * ierr_m(jj-1,:)'); % B_eval(:,:,jj) * (ddq_r(jj,:)'
        %we don't have the derivative error because we want that our
        %the derivative part of the controller acts as a pure damper

        %the torque for compensation of the robot inertia
        tau_comp(jj,:) = f*dq_msr(jj-1,:)' + g_eval(jj,:)';
        
        %this is the inertia matrix at the instant jj
        Beq = B(:,:,jj);
        
        %measured acceleration
        ddq_msr(jj,:) = (Beq)\( - tau_l(jj,:)' + tau_PID(jj,:)' + tau_comp(jj,:)');
        
        %measured speed
        dq_msr(jj,:) = dq_msr(jj-1,:) + ddq_msr(jj,:)*Ts;
        
        %measured position
        q_msr(jj,:) = q_msr(jj-1,:) + dq_msr(jj,:)*Ts;
        
        %integral error
        ierr_m(jj,:) = ierr_m(jj-1,:) + (q_r(jj,:) - q_msr(jj,:)) * Ts;
        
        %position error
        qerr(jj,:)   = q_r(jj,:) - q_msr(jj,:);

        %speed error
        dqerr(jj,:)  = dq_r(jj,:) - dq_msr(jj,:);

        %acceleration error
        ddqerr(jj,:) = ddq_r(jj,:) - ddq_msr(jj,:);
        
        %cost function J (idea: check penalty for instability)
        
        % delta_q = abs(q_msr(jj,:) - q_msr(jj-1,:));
        % max_delta_q = max(delta_q);
        
        if(abs(max(qerr(jj,:),[],2)) >= toll_qerr)
            %J = 1000;

            % now we try to implement a weighted penalty
            J = 1000*(1 - ((1 - p)/(log10(1 + length(t)))*log10(1 + jj)));
            exit_flag = true;
            max(qerr(jj-15:jj+5,:),[],2)

        end

        
        waitbar(jj/length(t),wb);
        
    end
    
    figure(15)
    title("qerr over iterations 1st-joint")
    plot(qerr(:,1))
    hold on

    %cost function J (idea: define KPIs for its implementation)
    J = J + sum(max(abs(qerr))) + sum(max(abs(dqerr))) + sum(mean(abs(qerr))) + sum(mean(abs(dqerr)));

    
    %close the waitbar
    close(wb)
    

    if (exit_flag == 1)
          %disp J terms
          fprintf('\npenalty: yes\n')
    else
        fprintf('\npenalty: no\n')
    end
    
    fprintf('Function evaluation %.0f: final cost: %12.8f \n', idx_sim, J)
    fprintf('--------------------------------\n')
    idx_sim = idx_sim + 1;
    
end