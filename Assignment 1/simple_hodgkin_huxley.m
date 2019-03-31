function simple_hodgkin_huxley
    % Adapted code originally written by Dr Eric Cook 
    % Modified by Thomas Christinck

    % Load the data
    load('CookAssignemnt1UnknownCurrent.mat');

    % The delta t in ms (time between updates)
    dt = .25;

    % Relevant info from data
    vstep_size = size(vStep);
    tstop = vstep_size(1) * dt;
    nb_voltages = vstep_size(2);

    %  Model parameters
    ba = .19;
    bi = -0.11;
    v0a = -43; 
    v0i = -31; 
    gBar = 0.001; % mS
    
    % The instantaneous current should be zero at the reversal potential
    Er = 0; % mV

    % We can figure out tau directly from the data
    taua = 3; %msec
    taui = 14; %msec

    % Plot steady state activation
    figure(1);
    clf;
    
    subplot(5,2,1);
    v = -80 : 20;
    
    % Plotting activation that is a function of ba (the slope at midpoint) and 
    % v0a (the voltage at midpoint). Similarly, we'll plot the inactivation function,
    % which is a function of bi and v0i.
    plot(v, two_param_sigmoid(v, [ba v0a]), 'LineWidth', 1.5);
    hold on;
    plot(v, two_param_sigmoid(v, [bi v0i]), 'LineWidth', 1.5);
    hold off;
    xlabel('mV');
    legend('xa(t = inf)', 'xi(t = inf)', 'Location', 'E');
    title('Activation and Inactivation Functions');
    
    % Plot parameters
    subplot(5,4,4);
    set(gca, 'Visible', 'off');
    text(0,0,['ba = ' num2str(ba)]);
    text(0,1,['v0a = ' num2str(v0a) ' mV']);
    text(0,2,['bi = ' num2str(bi)])
    text(0,3,['v0i = ' num2str(v0i) ' mV']);
    text(2,0,['taua = ' num2str(taua) ' msec']);
    text(2,1,['taui = ' num2str(taui) ' msec']);
    text(2,2,['Er = ' num2str(Er) ' mV']);
    text(2,3,['gBar = ' num2str(gBar) ' mS']);
    ylim([0 6]);
    xlim([0 4]);    
    
    
    % Voltages (there are 5 provided)
    v = zeros(nb_voltages, vstep_size(1));
    for k = 1:nb_voltages
        v(k, :) = vStep(:, k);
    end

    % Similarly, we set up the activation variables, inactivation variables, and conductance
    xa = zeros(nb_voltages, vstep_size(1));
    xaInf = zeros(nb_voltages);
    xi = zeros(nb_voltages, vstep_size(1));
    xiInf = zeros(nb_voltages);
    g = zeros(nb_voltages, vstep_size(1));
    
    % The delta t in ms (time between updates)
    dt = .25;
    
    % Compute the exponential terms (with constant taua, taui)
    ea = 1 - exp(-dt/taua);
    ei = 1 - exp(-dt/taui);

    tt = 0;
    j = 1;
    
    while tt < tstop
        % Now, for every update we make we'll have to consider each of the five
        % voltage curves given to us.

        % Update activation variable xaInf, inactivation variable xiInf
        for k = 1:nb_voltages
            xaInf(k) = two_param_sigmoid(v(k, j), [ba v0a]);
            xiInf(k) = two_param_sigmoid(v(k, j), [bi v0i]);
        end

        if j > 1
            % Each timestep, we have to update our xa, xi as the voltage
            % changes
            for k = 1:nb_voltages
                xa(k, j) = xa(k, j-1) + (xaInf(k) - xa(k, j-1)) * ea;
                xi(k, j) = xi(k, j-1) + (xiInf(k) - xi(k, j-1)) * ei;
            end
        else
            % At the start, the number of channels open is the steady-state
            % number (same for inactivation)
            for k = 1:nb_voltages
                xa(k, j) = xaInf(k);
                xi(k, j) = xiInf(k);
            end
        end
        
        % Update conductance and current 
        for k = 1:nb_voltages
            g(k, j) = gBar * xa(k, j) * xi(k, j);
            i(k, j) = g(k, j) * (v(k, j) - Er);
        end
        
        % Update time and then repeat procedure for computing conductance & current
        t(j) = tt;
        tt = tt + dt;
        j = j + 1;
    end
    
    % Plot time course
    subplot(6,1,2);
    plot(t, v, 'LineWidth', 1.5);
    axis tight;
    ylabel('V step (mV)');
    
    subplot(6,1,3);
    plot(t,xa, 'LineWidth', 1.5);
    axis tight;
    ylabel('xa');

    subplot(6,1,4);
    plot(t,xi, 'LineWidth', 1.5);
    axis tight;
    ylabel('xi');
    
    subplot(6,1,5);
    plot(t,g, 'LineWidth', 1.5);
    axis tight;
    ylabel('g (mS)');
    ylim([0 inf]);
    
    % Want to plot both target currents and model currents
    subplot(6,1,6);
    colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880]];
    for k = 1:nb_voltages
        % Plot the kth model current
        plot(t,i(k, :), 'LineWidth', 1.5, 'Color', colors(k, :));
        hold on;
        % Plot the kth target current (with dashed line)
        plot(t,iUnknownCurrent(:, k), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', colors(k, :));
        hold on;
    end
    hold off;
    axis tight;
    xlabel('msec');
    ylabel('i (mA)');
    ylim([-inf inf]);
end



