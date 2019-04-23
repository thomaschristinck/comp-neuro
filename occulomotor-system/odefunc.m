function dtheta = odefunc(t, y, F, p)
% ODE representing equation 2, dy/dt = (F - k * y) / r. Here y is theta, eye position, The force value for each
% time t is estimated by finding the time index closest to time t (interp1).

F = interp1(p.t, F, t);
dtheta = (F' - p.k * y) ./ p.r;

return;