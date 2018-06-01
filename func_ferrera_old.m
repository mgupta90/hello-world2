function [Ts_new, Ta_new, eta_new] = func_ferrera(Ta,Ts,F)

    % Model parameters
    yr = 360*24*3600; % One year in seconds
    a = 6.373*10^6; % Earth radius (m)
    Tf = -2; % Freezing temperature
    a0 = 0.72; % open ocean co-albedo
    a2 = -0.078; % Open ocean co-albedo (zenith angle dependence)
    Ai = 0.4; % Ice-covered albedo
    Ca = 10^7; % column integrated heat capcity for atmosphere [J/m2/C]
    Co = 18.5*10^7; % ocean heat capacity (MLD = 45m) [J/m2/C]
    Ka = 4.0 * 10^6; % atmospheric diffusivity [m2/s]
    Ko = 300; % oceanic diffusivity [m2/s]
    Aout = 202; % W/m2
    Bout = 1.9; % W/m2/K
    Aup = 238; % surface heat flux constant [W/m2]
    Bup = 15; % surface heat flux sensitiivty [W/m2/C]

    S0 = 1360 + F; % solar constant [W/m2]

    s2 = -0.48;
    dz = 50; % Depth scale for overturning (not sure about this value)

    % Space parameters 
    phi_max = pi/2;
    nlat = 50;
    dphi = phi_max/nlat;
    phi = dphi/2:dphi:(phi_max-dphi/2);
    phi_deg = phi * 180/pi; 

    P2_func = @(x_var) (3*x_var.^2-1)/2;
    A0 = a0 + a2*P2_func(sin(phi));
    S = S0/4*(1 + s2*P2_func(sin(phi)) ); % solar flux [W/m2]
    a_func = @(T_var) (T_var>Tf) .* A0 + (T_var<=Tf) .* Ai;

    PHI_res = zeros(1,nlat);
    PHI_res(1:nlat/2) = sin(phi(1:nlat/2)*4);

    % Time parameters
    tmax = 10; % Max simulation time in years
    dt = yr/10000; % timestep (seconds)
    t_int = 1; % Time interval between plots (years)

    % Declaring arrays
    t_plot = zeros(1,tmax/t_int+1);
    Ta_new = zeros(1,nlat);
    Ta_bar = zeros(1,tmax/t_int+1);
    Ts_new = zeros(1,nlat);
    Ts_bar = zeros(1,tmax/t_int+1);
    eta_deg = zeros(1,tmax/t_int+1);

    % % Initial conditions 
    load('Ts_ferrera4_warm.mat')
    load('Ta_ferrera4_warm.mat')

    Ta_prev = Ta;
    Ts_prev = Ts;

    [~,phicrit_idx] = min(abs(phi-50*pi/180)); %index of closest value
    Td = max(Ts_prev(phicrit_idx), Tf);
 
    eta_prev = get_eta(Ts_prev,Tf,phi_deg);
    eta_new = eta_prev + 5; % just to make sure it does at least one loop

    % Numerical integration
    num_it = 1;
    nsteps = tmax*yr/dt; % Total number of timesteps
    
%     while abs(eta_prev-eta_new) > 10^-15

    for i=2:nsteps+1
        eta_prev = eta_new;

        dTadphi = [0 diff(Ta_prev)./diff(phi) 0];
        dTsdphi = [0 diff(Ts_prev)./diff(phi) 0];

        Ts_mean = [min(Ts_prev) (Ts_prev(1:end-1) + Ts_prev(2:end))/2 max(Ts_prev)];
        phi_mean = [0 (phi(1:end-1) + phi(2:end))/2 phi_max];
        PHI_res_mean = [0 (PHI_res(1:end-1) + PHI_res(2:end))/2 0];
    

        Ha = -2*pi/a^2*Ca*Ka * cos(phi_mean) .* dTadphi;
        Hfa = -1./(2*pi*cos(phi)) .* diff(Ha) ./ diff(phi_mean);

        Ho = 2*pi*a*Co*cos(phi_mean) .* ( PHI_res_mean.*(Ts_mean-Td)./dz - Ko/a*dTsdphi);
        Hfo = -1./(2*pi*cos(phi))/a^2 .* diff(Ho) ./ diff(phi_mean);

        Fout = Aout + Bout*Ta_prev;
        alpha = a_func(Ts_prev);
        Sol = S.*alpha;
        Fup = Aup + Bup*(Ts_prev - Ta_prev);

        Fa_net = Hfa + Fup - Fout;
        Fs_net = Hfo - Fup + Sol;

        % New temperature profiles
        Ta_new = Ta_prev + dt/Ca*Fa_net;
        Ta_prev = Ta_new;
        Ts_new = Ts_prev + dt/Co*Fs_net;

        % Don't let Ts (SST) fall below Tf
        Ts_new(Ts_new < Tf) = Tf;

        Ts_prev = Ts_new;
        Td = max(Ts_prev(phicrit_idx), Tf);

        % New ice line
         eta_new = get_eta(Ts_prev,Tf,phi_deg);

         num_it = num_it + 1;

    end
    
    num_it

end


function eta = get_eta(Ts,Tf,phi)

    % Finding ice line
    [~,idx] = min(abs(Ts-Tf));
    
    if idx ~= 1 && idx ~= length(Ts)
        
        % Making sure the index is on the equatorial side (warmer)
        if abs(Ts(idx+1) - Ts(idx)) < abs(Ts(idx) - Ts(idx-1))
            idx = idx - 1; 
        end
        
        eta = phi(idx) + (phi(idx+1)-phi(idx))/(Ts(idx+1)-Ts(idx))*(Ts(idx+1)-Tf);
        
    
    
    else
        eta = phi(idx);
    end
end

