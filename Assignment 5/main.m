% A function that calculates occulomotor force based on eye-rotation amplitude and saccade duration (in part 1),
% by assuming constant force and solving the equation F = kθ + r * dθ/dt
% In part two, greater complexity is added. Robinson argued that the brain can adjust pulse frequency and amplitude
% via a simple feedback loop that produces an error signal = desired eye position (θdes) – actual estimated eye 
% position (θest). The error signal drives a high frequency (in neural terms) burst generator (BG) whose output is
% proportional to dθ/dt, which drives a circuit that compensates for the plant characteristics. There are 2 integrators 
% that perform the same function. However, to calculate error for each saccade the one that is used to calculate θest 
% needs to be reset to zero after every saccade. The output of the burst generator loop produces a firing frequency (ff) 
% in motorneurons -> this ff is proportional to the muscle force (assignment designed by D. Guitton, MNI)

% Completed independently by Thomas Christinck

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The given positions and durations
positions = [10 20 40];
durations = [0.045 0.068 0.110];

%Stiffness and viscosity (k and r) parameters from Sylvester and Cullen (1999)
%Parameters----------------------------------------------------------------
p.r  = 0.42;        %g.sec/deg: viscous element
p.k  = 4.2;        %g/deg:     spring stiffness
p.dt = 0.005;       %sec:       sampling time
p.t  = 0:p.dt:0.3;  %time:      time vector
sac_dur = 110e-3;   %sec:       saccade duration
sac_mag = 40;       %deg:       saccade magnitude
%--------------------------------------------------------------------------

% For each position/duration, find the force during the movement as well as after. Then plot theta as well as the force
for i = 1:3
	F = p.k * positions(i) / (1 - exp((-p.k * durations(i)) / p.r));
	fprintf("\nForce during saccade for position %d and duration %.3f is %.2f \n", positions(i), durations(i), F);
	F_post = p.k * positions(i);
	fprintf("Force after saccade for position %d and duration %.3f is %.2f \n", positions(i), durations(i), F_post);
	% We multiply by 1000 to plot in ms
	theta = @(t) (F / p.k * (1 - exp(- p.k * (t) / p.r))).*(t < durations(i)) + (positions(i)).*(t >= durations(i));
	figure(1)
	fplot(theta, [0 0.150], 'Linewidth', 1.5); hold on;
	F = @(t) (F).*(t < durations(i)) + F_post.*(t >= durations(i));
	figure(2)
	fplot(F, [-0.0005 0.140], 'Linewidth', 1.5); hold on;

end

% Set up plots for each of the saccades; both for amplitude and force
figure(1)
title('Eye Rotation Amplitude as a Function of Time')
xlabel('Time (s)')
ylabel('Eye Rotation (degrees)')
legend('10^{\circ} (45ms)', '20^{\circ} (68ms)', '40^{\circ} (110ms)', 'Location', 'nw')

figure(2)
title('Force as a Function of Time')
xlabel('Time (s)')
ylabel('Force (g)')
legend('10^{\circ} (45ms)', '20^{\circ} (68ms)', '40^{\circ} (110ms)', 'Location', 'ne')
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the control system to find the estimated position, BG Output (proportional to velocity), and
% force as a function of time
[theta_est, velocity, force] = runControlSys(sac_mag, sac_dur, p);
figure(3)
plot(p.t, theta_est, p.t, velocity, p.t, force, 'Linewidth', 1.5)
legend('Eye position (deg)', 'Velocity (deg/s)', 'Force (g)')
xlabel('Time (s)')
title('Estimated Eye Position, Burst Generator Output, and Force for a 40^{\circ}, 110 ms Saccade ')
ylabel ('Magnitude')

% Solve the ODE from equation 2
[t,theta_out] = ode45(@(t,theta_out) odefunc(t,theta_out,force,p), p.t, 0);

% Plot theta from part one (1st order ODE only), as well as theta estimated from the neural circuit and the
% actual theta from the eye plant
figure(4)
fplot(theta, [-0.0005 0.140], 'Linewidth', 1.2); hold on;
plot(p.t, theta_est, 'Linewidth', 1.2);
plot(p.t, theta_out, 'Linewidth', 1.2);
xlim([0 0.14])
title('Neural Circuit-Estimated and Plant-Output Eye Rotation for 40^{\circ}, 110 ms Saccade')
xlabel('Time (s)')
ylabel('Position (degrees)')
legend('Constant force eye position', 'Feedback loop estimated eye position \theta_{est}', 'Plant output eye position \theta_{out}', 'Interpreter', 'tex', 'Location', 'se')

% Now, we assume that the integrator output in the loop (Box 1) is deficient and that its gain, k_lesion, 
% is half of what it should be to compensate for the plant’s (Box 2) k. 
p.k = p.k / 2;
[theta_est,veldelta,force] = runControlSys(sac_mag, sac_dur, p);   
p.k = p.k * 2;
[t,theta_out_lesion] = ode45(@(t,theta_out) odefunc(t,theta_out,force,p), p.t, 0);

%Show the resultant eye trajectories, θ_out for a 40° saccade. (Hint: First calculate the loop output force 
% temporal profile for klesion and then calculate θout using the plant equation and the force profile.) 
figure(5)
plot(p.t, theta_out_lesion, p.t, theta_out, 'Linewidth', 1.5);
title('40^{\circ}, 110 ms Saccade Eye Trajectory for Deficient and Non-Deficient Neural Circuits')
xlabel('Time (s)')
ylabel('Position (degrees)')
legend('Deficient Control Circuit with k_{lesion} = 0.5 * k', 'Non-Deficient Control Circuit', 'Interpreter', 'tex', 'Location', 'se');




