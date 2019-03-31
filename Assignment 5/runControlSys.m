function [theta_est, velocity, force]=run_ControlSys(sac_mag, sac_dur, p)
% Function that computes the estimated eye position signal (θest), the BG output and the Force (ff) 
% versus time for a 40 degree saccade; i.e. θdes = 40 deg. (Set the BG output such that eye reaches 
% saccade-end in 110ms.)

% Params:
% magnitude - the saccade magnitude in degrees
% duration - the saccade duration
% time - a time vector for integration
% dt - time step for integration
% r, k - parameters for viscosity and stiffness

theta_des = sac_mag;           	% saccade magnitude (deg)
theta_est = zeros(1,length(p.t));  % eye position (deg)
velocity = zeros(1,length(p.t));   % eye velocity (deg/s): 
force  = zeros(1,length(p.t));    	% force (g)
BG = sac_mag / sac_dur;         	% burst generator output (deg/s)

% Now, for each time we compute the error - if error is non-zero then the burst generator is constant, if
% error is zero then burst generator is silenced
for n = 2:length(p.t)
    error = theta_des - theta_est(n-1);
    if error > 0
        velocity(n) = BG;
    else
        velocity(n) = 0;
    end
    % Compute the integral under velocity curve
    theta_est(n) = trapz(1:n, velocity(1:n) .* p.dt);
    % Equation in box 1 gives us the force
    force(n) = p.k * theta_est(n) +  p.r * velocity(n);
end

return;
