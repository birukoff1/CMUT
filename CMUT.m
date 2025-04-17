classdef CMUT

    properties

        R_pl % Plate radius 

        H_pl % Thickness of plate
        D_w % Pillars width
        H_ins_l % Thickness of insulator layer
        H_gap % Thickness of gap


        C % Capacitance
        Q % Charge
        
        D_elec % Averaged electric displacement field
        F % Electrostatic force
        D_rigid % Bending stiffness

        Disp_max % Maximum displacement in a static mode
        Disp_static % Static displacement profile
        V_collapse % Collapse voltage

        f_1 % Frequency of the first harmonic
        f_1_shifted % Frequency of the first harmonic (dc shift)
        
        Eigenfrequencies % Eigenfrequencies
        Eigenmodes % Displacement field of the eigenmodes
        Amplitude_Response % Amplitude response

        v_n % Normal velocity

        Frequency_Response % Frequency response
        Radiation_pattern % Radiation pattern
    end


    properties (Hidden = true)

        R % Effective radius
        S % Effective area
        
        r
        theta

        z
        phi

        t

        V_dc % Applied DC voltage
        H_gap_initial

        rho_medium = 1000; % Density of the medium

    end


    properties (Constant, Hidden = true)

        Epsilon_0 = 8.85e-12;
        Epsilon_si = 4.5;
        Epsilon_si_ni = 9.7;
        Epsilon_vac = 1.0058;


        E_si = 160e9;
        rho_si = 2320;
        mu_si = 0.22;

        E_si_ni = 250e9;
        rho_si_ni = 3100;
        mu_si_ni = 0.23;

        E_al = 70e9;
        rho_al = 2700;
        mu_al = 0.35;

        N_r = 101; % Number of points in a radial direction
        N_theta = 180; % Number of points in an angular direction

        N_phi = 181; % Number of angles for radiation pattern

        N_t = 100; % Number of time points

    end



    methods

        function obj = CMUT(R_pl, H_pl, D_w, H_ins_l, H_gap)
            
            % Assigning values to geometrical parameters
            obj.R_pl = R_pl*1e-6;
    
            obj.H_pl = H_pl*1e-6;
            obj.D_w = D_w*1e-6;
            obj.H_ins_l = H_ins_l*1e-6;
            obj.H_gap = H_gap*1e-6;
            obj.H_gap_initial = H_gap*1e-6;
            
            obj.R = obj.R_pl-obj.D_w;
            obj.S = pi*obj.R^2;

            
            % Calculation of constant mechanical parameters
            obj.D_rigid = (obj.E_si*obj.H_pl^3)/(12*(1-obj.mu_si^2));


            % Calculation of the collapse voltage
            k = (16*pi*obj.E_si*obj.H_pl^3)/(3*(1-obj.mu_si^2)*obj.R^2); % Stiffness coefficient
            H_eff = obj.H_ins_l/obj.Epsilon_si_ni+obj.H_gap_initial/obj.Epsilon_vac;
            obj.V_collapse = 1.2*sqrt((8*k*H_eff^3)/(27*obj.Epsilon_0*obj.S));


            % Specifying the spatial and temporal variables
            obj.r = linspace(-obj.R, obj.R, obj.N_r)';
            obj.theta = linspace(0,359, obj.N_theta)*pi/180;
            obj.t = linspace(0, 1e8, obj.N_t)';
            
        end


        function obj = Find_Lumped_Parameters(obj, V_dc)
        
            % Calculation of constant electrical parameters (the pillars are taken into consideration)
            obj.C = obj.Epsilon_0*obj.S*(1/(obj.H_ins_l/obj.Epsilon_si_ni+obj.H_gap_initial/obj.Epsilon_vac)) + obj.Epsilon_0*(pi*obj.R_pl^2-obj.S)*(1/((obj.H_ins_l+obj.H_gap_initial)/obj.Epsilon_si_ni));
            
            % Calculation of variable electrical parameters 
            obj.Q = obj.C*V_dc;
            obj.V_dc = V_dc;

            if obj.V_dc > 0.9*obj.V_collapse
                disp('Danger! DC voltage is higher than 90% of the collapse voltage')
            end

            
            % Calculation of static deformation
            Disp_max_new = 1;
            obj.Disp_max = 10;

            Disp_new = zeros(obj.N_r, 1);
            F_elec = zeros(obj.N_r, 1);

            N_r_new = 1e3;
            r_new = linspace(-obj.R, obj.R, N_r_new)';
            dr = (2*obj.R)/(N_r_new-1);
            F_mec = zeros(N_r_new, N_r_new);

            for i=3:N_r_new-2
                if (i~=1) && (i~=2)
                    F_mec(i,i-2) = 1/dr^4 - 2/r_new(i)*1/(2*dr^3) - 1/r_new(i)^2*0/dr^2 + 1/r_new(i)^3*0/(2*dr);
                end
                if (i~=1)
                    F_mec(i,i-1) = -4/dr^4 + 2/r_new(i)*2/(2*dr^3) - 1/r_new(i)^2*1/dr^2 - 1/r_new(i)^3*1/(2*dr);
                end
                F_mec(i,i) = 6/dr^4 + 2/r_new(i)*0/(2*dr^3) + 1/r_new(i)^2*2/dr^2 + 1/r_new(i)^3*0/(2*dr);
                if (i~=N_r_new)
                    F_mec(i,i+1) = -4/dr^4 - 2/r_new(i)*2/(2*dr^3) - 1/r_new(i)^2*1/dr^2 + 1/r_new(i)^3*1/(2*dr);
                end
                if (i~=N_r_new) && (i~=N_r_new-1)
                    F_mec(i,i+2) = 1/dr^4 + 2/r_new(i)*1/(2*dr^3) - 1/r_new(i)^2*0/dr^2 + 1/r_new(i)^3*0/(2*dr);
                end
            end
            F_mec = F_mec(:,3:end-2);

            n = 0;
            while ( abs(Disp_max_new-obj.Disp_max)/Disp_max_new > 1e-3 || n < 2 )
                    
                obj.Disp_max = Disp_max_new;
                Disp = Disp_new;

                F_elec = 1/2*obj.Epsilon_0./((obj.H_gap_initial-Disp)/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*V_dc^2;

                if n == 0
                    F_elec = 2*pi*trapz(obj.r(ceil(obj.N_r/2):end), F_elec(ceil(obj.N_r/2):end).*obj.r(ceil(obj.N_r/2):end));
                    obj.Disp_max = (F_elec/obj.S*obj.R^4)/(64*obj.D_rigid);
                    Disp_new = obj.Disp_max/obj.R^4.*(obj.R^2-obj.r.^2).^2;
                else

                    F_elec = interp1(obj.r, F_elec, r_new);
                    Disp_new = linsolve(F_mec, F_elec/obj.D_rigid);
                    Disp_new = [0; Disp_new(1)/2; Disp_new; Disp_new(end)/2; 0];
                    Disp_new = interp1(r_new, Disp_new, obj.r);
                end
                Disp_max_new = max(Disp_new);

                n=n+1;
                if Disp_max_new > obj.H_gap
                    disp('Short-circuit! Tht DC voltage must be decreased')
                    error
                end
            end
            if n==1
                obj.F = F_elec;
            else
                obj.F = 2*pi*trapz(r_new(ceil(N_r_new/2):end), F_elec(ceil(N_r_new/2):end).*r_new(ceil(N_r_new/2):end));
            end
            obj.Disp_max = Disp_max_new;
            obj.Disp_static = Disp_new;
            

            % Calculation of eigenfrequencies
            H_eff = obj.H_ins_l/obj.Epsilon_si_ni+obj.H_gap_initial/obj.Epsilon_vac;
            k = (16*pi*obj.E_si*obj.H_pl^3)/(3*(1-obj.mu_si^2)*obj.R^2);
            k_shifted = k - 1*obj.Epsilon_si*obj.S*obj.V_dc^2/H_eff; % Coefficient 1
            
            obj.f_1 = 1/(2*pi)*sqrt(k/(0.613*obj.S*obj.H_pl*obj.rho_si)); % Coefficient from mass to effective mass (0.613)
            obj.f_1_shifted = 1/(2*pi)*sqrt(k_shifted/(0.613*obj.S*obj.H_pl*obj.rho_si));
            % obj.f_1 = 5.1/(pi*obj.R^2)*sqrt(obj.D_rigid/(obj.rho_si_ni*obj.H_pl));
        end


        function obj = Find_Eigenvalues(obj)

            m = obj.H_pl*obj.rho_si; % Mass per area kg/m^2
            
            Beta = [3.19622   4.61090   5.90568   7.14353;      6.30644   7.79927   9.19688  10.53667;      9.43950  10.95807  12.40222  13.79506; ...
                    12.57713  14.10863  15.57949  17.00529;     15.71644  17.25573  18.74396  20.19231;     18.85655  20.40104  21.90149  23.36628]'; % Solution of the eigenproblem for a circular plate

            N_c = length(Beta(1,:)); % Number of circles
            N_d = length(Beta(:,1)); % Number of diameters

            obj.Eigenfrequencies = zeros(N_d, N_c);
            obj.Eigenmodes = zeros(length(obj.r), length(obj.theta), N_d, N_c);

            A_kn = ones(N_d, N_c); % Coefficients of contribution of the eigenmodes
            for i=1:N_d
                for j=1:N_c
                    
                    obj.Eigenfrequencies(i,j) = 1/(2*pi)*sqrt((Beta(i,j)/obj.R)^4*obj.D_rigid/m); % Beta^2 is proportional to eigenvalue
                    obj.Eigenfrequencies(i,j) = obj.Eigenfrequencies(i,j)*obj.f_1_shifted/obj.f_1; % Adjustment in accordance to the calculated f_1_shifted

                    obj.Eigenmodes(:,:,i,j) = (besseli(i-1,Beta(i,j))*besselj(i-1,Beta(i,j)*obj.r/obj.R) - ...
                        besselj(i-1,Beta(i,j))*besseli(i-1,Beta(i,j)*obj.r/obj.R))*cos((i-1)*obj.theta);
                    
                    I = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), m*obj.Eigenmodes(ceil(obj.N_r/2):end,:,i,j).^2.*obj.r(ceil(obj.N_r/2):end), 1));
                    A_kn(i,j) = sqrt(1/I); % Normalization coefficient

                    obj.Eigenmodes(:,:,i,j) = obj.Eigenmodes(:,:,i,j)*A_kn(i,j); % Setting an amplitude
                end
            end

            % Finding coefficients of the eigenmodes
            q_dc = 1/2*obj.Epsilon_0./(obj.H_gap/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*obj.V_dc^2;
            Factor_dc = zeros(length(obj.Eigenfrequencies(1,:)),1);
            Mode_response = zeros(length(obj.r), length(Factor_dc) );
            for i=1:length(Factor_dc)

                Factor_dc(i) = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), obj.Eigenmodes(ceil(obj.N_r/2):end,:,1,i)*q_dc.*obj.r(ceil(obj.N_r/2):end), 1));
                Mode_response(:,i) = obj.Eigenmodes(:,1,1,i)/(2*pi*obj.Eigenfrequencies(1,i))^2*Factor_dc(i);
            end
            Response = sum(Mode_response(ceil(obj.N_r/2),:));
            obj.Eigenmodes = obj.Eigenmodes*sqrt(obj.Disp_max/Response); % Coefficients are set so the static response is equal to calculated directly
        end


        function obj = Plot_Eigenmodes(obj, j, k, dimension)

            if (j > length(obj.Eigenfrequencies(1,:)))
                disp(['Only ', num2str(length(obj.Eigenfrequencies(1,:))-1), ' node circles are available'])
                error
            elseif (j < 1)
                disp('Minimum number of node circles is 1')
                error
            elseif (k > length(obj.Eigenfrequencies(:,1))-1)
                disp(['Only ', num2str(length(obj.Eigenfrequencies(:,1))-1), ' node diameters are available'])
                error
            elseif (k < 0)
                disp('Minimum number of node diameters is 0')
                error
            end


            if strcmp(dimension,"1D") == 1
                plot(1e6*obj.r, 1e6*1/5*obj.R*obj.Eigenmodes(:,1,k+1,j)) % 1/5 is a scaling coefficient
                title(['Eigenmode on a frequency = ',num2str(1e-6*obj.Eigenfrequencies(k+1,j+1)), ' MHz'], ['Amount of node diameters = ',num2str(k), ', Amount of node circles = ',num2str(j)])
                xlim([-1e6*obj.R; 1e6*obj.R])
                xlabel('Radial coordinate, um')
                ylabel('Amplitude, um')
                f = gca;
                %exportgraphics(f,'Eigenmode_1D.png','Resolution',300)

            elseif strcmp(dimension,"2D") == 1
                [theta_mesh,r_mesh] = meshgrid(obj.theta, obj.r);
                [Theta_polar,R_polar,Plate_disp_polar] = pol2cart(theta_mesh, r_mesh, obj.Eigenmodes(:,:,:,j));
                mesh(1e6*Theta_polar,1e6*R_polar, 1e6*1/4*obj.R*Plate_disp_polar(:,:,k+1)/max(Plate_disp_polar(:,:,k+1), [], 'all'));
                title(['Eigenmode on a frequency = ',num2str(1e-6*obj.Eigenfrequencies(k+1,j)), ' MHz'], ['Amount of node diameters = ',num2str(k), ', Amount of node circles = ',num2str(j)])
                xlabel('Radial coordinate, um')
                zlabel('Amplitude, um')
                colormap("cool")
                shading interp
                axis equal
                f = gca;
                %exportgraphics(f,'Eigenmode_2D.png','Resolution',300)

            else
                disp('"1D" or "2D" plots are supported')
            end
        end


        function obj = Find_Response(obj, V_dc, V_ac, f)
            
            V_dc = abs(V_dc);
            V_ac = abs(V_ac);

            f = 1e6*f; % From MHz to Hz
            if V_dc > 0.9*obj.V_collapse
                disp('Danger! DC voltage is higher than 90% of the collapse voltage')
            end

            if V_dc~=obj.V_dc
                disp('New DC voltage value')
                disp('Update of lumped parameters and eigenfrequencies')
                obj = Find_Lumped_Parameters(obj, V_dc, load('input_Mean.txt'));
                obj = Find_Eigenvalues(obj);
            end


            % Static case, when V_ac = 0
            if (V_dc ~= 0) && (V_ac == 0)

                Disp = zeros(obj.N_r, 1);
                F_dc = 1/2*obj.Epsilon_0*obj.S./((obj.H_gap-Disp)/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*V_dc^2;

                q_dc = zeros(obj.N_r, obj.N_theta);
                for i=1:obj.N_theta
                    q_dc(:,i) = F_dc/obj.S;
                end

                Factor_dc = zeros(length(obj.Eigenfrequencies(1,:)),1); % Participation factors
                Mode_response = zeros(length(obj.r), length(Factor_dc) ); % Response of every mode to external excitation f(r,n)
                for i=1:length(Factor_dc)
    
                    Factor_dc(i) = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), obj.Eigenmodes(ceil(obj.N_r/2):end,:,1,i).*q_dc(ceil(obj.N_r/2):end,:).*obj.r(ceil(obj.N_r/2):end), 1)); % Calculation of the participation factors
                    omega_i = 2*pi*obj.Eigenfrequencies(1,i);
                    Mode_response(:,i) = obj.Eigenmodes(:,1,1,i)/omega_i^2*Factor_dc(i); % Calculation of the solution to ODE
                end

                obj.Amplitude_Response.static.amp = zeros( length(obj.r), length(obj.theta) );
                for i=1:length(obj.r)
                    obj.Amplitude_Response.static.amp(i,:) = sum(Mode_response(i,:));
                end


            % Dynamic case, when V_ac = A*sin(wt)
            elseif (V_dc ~= 0) && (V_ac ~= 0) && (length(f) == 1)
            
                Disp = zeros(obj.N_r, 1);
                F_dc = 1/2*obj.Epsilon_0*obj.S./((obj.H_gap-Disp)/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*V_dc^2;
                F_ac = 1/2*obj.Epsilon_0*obj.S./((obj.H_gap-Disp)/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*2*V_ac*V_dc;
                q_dc = zeros(obj.N_r, obj.N_theta); % Constant force, distributed along the rings
                q_ac = zeros(obj.N_r, obj.N_theta); % Alternative force, distributed along the rings
                for i=1:obj.N_theta
                    q_dc(:,i) = F_dc/obj.S;
                    q_ac(:,i) = F_ac/obj.S;
                end
    
                Factor_dc = zeros(length(obj.Eigenfrequencies(1,:)),1); % Participation factors
                Factor_ac = zeros(length(obj.Eigenfrequencies(1,:)),1);
    
                Mode_response_static = zeros(length(obj.r), length(obj.Eigenfrequencies(1,:)) ); % Static response of every mode to external excitation f(r,n)
                Mode_response_dynamic = zeros(length(obj.r), length(obj.Eigenfrequencies(1,:)) ); % Dynamic response of every mode to external excitation f(r,n)
                Mode_response = zeros(length(obj.r), length(obj.Eigenfrequencies(1,:)) ); % Full response of every mode to external excitation f(r,n)
    
                for i=1:length(Factor_dc)
                    Factor_dc(i) = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), obj.Eigenmodes(ceil(obj.N_r/2):end,:,1,i).*q_dc(ceil(obj.N_r/2):end,:).*obj.r(ceil(obj.N_r/2):end), 1)); % Calculation of the participation factors
                    Factor_ac(i) = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), obj.Eigenmodes(ceil(obj.N_r/2):end,:,1,i).*q_ac(ceil(obj.N_r/2):end,:).*obj.r(ceil(obj.N_r/2):end), 1));
                    
                    T1 = 1/obj.Eigenfrequencies(1,i);
                    T2 = 1/f;
                    omega_i = 2*pi*obj.Eigenfrequencies(1,i);
                    omega_exc = 2*pi*f;
                    order = 1e6/min(T1, T2);
                    T_new = lcm(round(order*T1), round(order*T2))/order;
                    t_new = linspace(0, T_new, 1e4);
                    w_new = abs(Factor_ac(i)/(omega_exc^2-omega_i^2)) * (cos(omega_exc*t_new) - cos(omega_i*t_new));
    
                    Mode_response_static(:,i) = obj.Eigenmodes(:,1,1,i)*Factor_dc(i)/omega_i^2; % DC Responses
                    Mode_response_dynamic(:,i) = obj.Eigenmodes(:,1,1,i)*sqrt(1/T_new*trapz(t_new,w_new.^2)); % AC Responses
                    Mode_response(:,i) = Mode_response_static(:,i) + Mode_response_dynamic(:,i); % DC + AC Responses
                end
                
                obj.Amplitude_Response.dynamic.factors = max(Mode_response_dynamic);
                obj.Amplitude_Response.dynamic.f = f;
                obj.Amplitude_Response.dynamic.amp = zeros( length(obj.r), length(obj.theta), 3 );
                for i=1:length(obj.r)
                    obj.Amplitude_Response.dynamic.amp(i,:,1) = sum(Mode_response_static(i,:));
                    obj.Amplitude_Response.dynamic.amp(i,:,2) = sum(Mode_response_dynamic(i,:));
                    obj.Amplitude_Response.dynamic.amp(i,:,3) = sum(Mode_response(i,:));
                end


            % Sweep case, when V_ac = A*sin([w_1 ... w_n]*t)
            elseif (V_dc ~= 0) && (V_ac ~= 0) && (length(f) ~= 1)

                Disp = zeros(obj.N_r, 1);
                F_dc = 1/2*obj.Epsilon_0*obj.S./((obj.H_gap-Disp)/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*V_dc^2;
                F_ac = 1/2*obj.Epsilon_0*obj.S./((obj.H_gap-Disp)/obj.Epsilon_vac+obj.H_ins_l/obj.Epsilon_si_ni).^2*2*V_ac*V_dc;
                q_dc = zeros(obj.N_r, obj.N_theta);
                q_ac = zeros(obj.N_r, obj.N_theta);
                for i=1:obj.N_theta
                    q_dc(:,i) = F_dc/obj.S;
                    q_ac(:,i) = F_ac/obj.S;
                end
   
                Factor_dc = zeros(length(obj.Eigenfrequencies(1,:)),1); % Participation factors
                Factor_ac = zeros(length(obj.Eigenfrequencies(1,:)),1);

                Mode_response_static = zeros(length(obj.r), length(obj.Eigenfrequencies(1,:)) ); % Static response of every mode to external excitation f(r,n)
                Mode_response_dynamic = zeros(length(obj.r), length(obj.Eigenfrequencies(1,:)) ); % Dynamic response of every mode to external excitation f(r,n)
                Mode_response = zeros(length(obj.r), length(obj.Eigenfrequencies(1,:)) ); % Full response of every mode to external excitation f(r,n)

                obj.Amplitude_Response.sweep.f = f;
                obj.Amplitude_Response.sweep.amp = zeros( length(obj.r), length(obj.theta), length(f), 3);

                for j=1:length(f)
                    for i=1:length(Factor_dc)
        
                    Factor_dc(i) = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), obj.Eigenmodes(ceil(obj.N_r/2):end,:,1,i).*q_dc(ceil(obj.N_r/2):end,:).*obj.r(ceil(obj.N_r/2):end), 1)); % Calculation of the participation factors
                    Factor_ac(i) = trapz(obj.theta, trapz(obj.r(ceil(obj.N_r/2):end), obj.Eigenmodes(ceil(obj.N_r/2):end,:,1,i).*q_ac(ceil(obj.N_r/2):end,:).*obj.r(ceil(obj.N_r/2):end), 1));

                    T1 = 1/obj.Eigenfrequencies(1,i);
                    T2 = 1/f(j);
                    omega_i = 2*pi*obj.Eigenfrequencies(1,i);
                    omega_exc = 2*pi*f(j);
                    order = 1e6/min(T1, T2);
                    T_new = lcm(round(order*T1), round(order*T2))/order;
                    t_new = linspace(0, T_new, 1e4);
                    w_new = Factor_ac(i)/(omega_exc^2-omega_i^2) * (cos(omega_exc*t_new) - cos(omega_i*t_new));

                    Mode_response_static(:,i) = obj.Eigenmodes(:,1,1,i)*Factor_dc(i)/omega_i^2; % DC Responses
                    Mode_response_dynamic(:,i) = obj.Eigenmodes(:,1,1,i)*sqrt(1/T_new*trapz(t_new,w_new.^2)); % AC Responses 
                    Mode_response(:,i) = Mode_response_static(:,i) + Mode_response_dynamic(:,i); % DC + AC Responses
                    end

                    for i=1:length(obj.r)
                        obj.Amplitude_Response.sweep.amp(i,:,j, 1) = sum(Mode_response_static(i,:));
                        obj.Amplitude_Response.sweep.amp(i,:,j, 2) = sum(Mode_response_dynamic(i,:));
                        obj.Amplitude_Response.sweep.amp(i,:,j, 3) = sum(Mode_response(i,:));
                    end  
                end

            else
                obj.Amplitude_Response.dynamic.factors(:) = 0;
                obj.Amplitude_Response.dynamic.amp(:,:) = 0;
            end
        end

       

        function obj = Plot_Response(obj, type, dimension)
            

            if strcmp(type,"Static") == 1
                
                if strcmp(dimension,"1D") == 1
                    figure
                    plot(1e6*obj.r,-1e6*obj.Amplitude_Response.static.amp(:,1), 1e6*obj.r,-1e6*obj.Disp_max/obj.R^4.*(obj.R^2-obj.r.^2).^2, '--')
                    title('Response to the static load')
                    xlabel('Radial coordinate, um')
                    ylabel('Amplitude, um')
                    xlim([1e6*obj.r(1); 1e6*obj.r(end)])
                    ylim([-1.05*1e6*max(obj.Disp_max, max(obj.Amplitude_Response.static.amp(:,1))); 0.01*1e6*obj.Disp_max])
                    legend('Mode superposition', 'Direct calculation', 'Location','southeast')
                    f = gca;
                    %exportgraphics(f,'StaticResponse_1D.png','Resolution',300)

                elseif strcmp(dimension,"2D") == 1

                    [theta_mesh,r_mesh] = meshgrid(1e6*obj.theta, 1e6*obj.r);
                    
                    figure
                    [Theta_polar,R_polar,Response_polar] = pol2cart(theta_mesh, r_mesh, 1e6*obj.Amplitude_Response.static.amp);
                    mesh(Theta_polar,R_polar, 1/4*obj.R*1e6*Response_polar/max(Response_polar, [], "all"))
                    title('Response to the static load')
                    xlabel('Radial coordinate, um')
                    zlabel('Amplitude, um')
                    colormap("cool")
                    shading interp
                    axis equal
                    f = gca;
                    %exportgraphics(f,'StaticResponse_2D.png','Resolution',300)
                else
                    disp('"1D" or "2D" plots are supported')
                end


            elseif strcmp(type,"Dynamic") == 1

                if strcmp(dimension,"1D") == 1

                    tiledlayout(2,1)

                    nexttile
                    plot(1e-6*obj.Eigenfrequencies(1,:), 1e6*obj.Amplitude_Response.dynamic.factors, '*')
                    title('Contribution of the modes')
                    xlabel('Eigenfrequency, MHz')
                    ylabel('Participation factor, -')
                    zlabel('Amplitude, um')
                    
                    nexttile
                    plot(1e6*obj.r,-1e6*obj.Amplitude_Response.dynamic.amp(:,1,3), 1e6*obj.r,-1e6*obj.Amplitude_Response.dynamic.amp(:,1,1), '--')
                    title('Response to the dynamic load')
                    xlabel('Radial coordinate, um')
                    ylabel('Amplitude, um')
                    legend('Full response', 'Static deflection')
                    ylim([-1.05*1e6*max(max(obj.Amplitude_Response.dynamic.amp(:,1,3)), max(obj.Amplitude_Response.dynamic.amp(:,1,1))); ...
                        0.01*1e6*max(max(obj.Amplitude_Response.dynamic.amp(:,1,3)), max(obj.Amplitude_Response.dynamic.amp(:,1,1)))])
                    xlim([1e6*obj.r(1); 1e6*obj.r(end)])
                    f = gca;
                    %exportgraphics(f,'DynamicResponse_1D.png','Resolution',300)

                elseif strcmp(dimension,"2D") == 1

                    [theta_mesh,r_mesh] = meshgrid(1e6*obj.theta, 1e6*obj.r);
                    
                    tiledlayout(2,1)

                    nexttile
                    plot(1e-6*obj.Eigenfrequencies(1,:), 1e6*obj.Amplitude_Response.dynamic.factors, '*')
                    title('Contribution of the modes')
                    xlabel('Eigenfrequency, MHz')
                    ylabel('Participation factor, -')
                    zlabel('Amplitude, um')
                    
                    nexttile
                    [Theta_polar,R_polar,Response_polar] = pol2cart(theta_mesh, r_mesh, 1e6*obj.Amplitude_Response.dynamic.amp(:,:,3));
                    mesh(Theta_polar,R_polar, 1/4*obj.R*1e6*Response_polar/max(Response_polar, [], "all"))
                    title('Response to the dynamic load')
                    xlabel('Radial coordinate, um')
                    zlabel('Amplitude, um')
                    colormap("cool")
                    shading interp
                    axis equal
                    f = gca;
                    %exportgraphics(f,'DynamicResponse_2D.png','Resolution',300)
                else
                    disp('"1D" or "2D" plots are supported')
                end


            elseif strcmp(type,"Sweep") == 1
                
                Response = zeros(length(obj.Amplitude_Response.sweep.f), 1);
                for i=1:length(obj.Amplitude_Response.sweep.f)
                    if strcmp(dimension,"Static") == 1
                        Response(i) = obj.Amplitude_Response.sweep.amp(ceil(obj.N_r/2),1,i, 1);

                    elseif strcmp(dimension,"Dynamic") == 1
                        Response(i) = obj.Amplitude_Response.sweep.amp(ceil(obj.N_r/2),1,i, 2);

                    elseif strcmp(dimension,"Full") == 1
                        Response(i) = obj.Amplitude_Response.sweep.amp(ceil(obj.N_r/2),1,i, 3);
                    else
                        disp('"Static", "Dynamic" or "Full" sweep response types are supported')
                    end

                end
                figure
                plot(1e-6*obj.Amplitude_Response.sweep.f, log10( 1e6*abs(Response)), '--' )
                title('Amplitude response')
                xlabel('External frequency, MHz')
                ylabel('Amplitude, log_{10}(um)')
                f = gca;
                %exportgraphics(f,'SweepResponse_1D.png','Resolution',300)

            else
                disp('"Static", "Dynamic" or "Sweep" types are supported')
            end
        end



        function obj = Find_FR(obj, type, Nx, Ny, phase, d, medium, r0, direction)
            
            r0 = r0*1e-6;
            % Setting parameters for the medium
            
            if strcmp(medium,"air") == 1
                c = 343;
                p_ref = 20e-6;
                obj.rho_medium = 1.2754;
            elseif strcmp(medium,"water") == 1
                c = 1500;
                p_ref = 1e-6;
                obj.rho_medium = 997;
            end


            % Assigning normal veocity to a grid of the array
            
            d = d*1e-6;
            if d<2*obj.R_pl
                disp(['Distance must be greater than a diameter of one cell: ', num2str(2*obj.R_pl*1e6), ' um'])
                error
            end

            x = (d*(Nx-1)+2*obj.R); % Length of the array in x direction
            y = (d*(Ny-1)+2*obj.R); % Length of the array in y direction

            Scale = 1;
            dr = obj.r(ceil(obj.N_r/2)+1)*Scale; % Grid step

            N_x = round(x/dr);
            if mod(N_x,2)==0
                N_x = N_x+1; % Odd number of points in x direction
            end

            N_y = round(y/dr);
            if mod(N_y,2)==0
                N_y = N_y+1; % Odd number of points in y direction
            end
            
            x = dr*(-floor(N_x/2):floor(N_x/2)); % Array of points in x direction
            y = dr*(-floor(N_y/2):floor(N_y/2)); % Array of points in y direction
          

            if r0 < max(x(end),y(end))^2 / (c/obj.Amplitude_Response.dynamic.f)
                disp('The distance from the array is too short. Set bigger distance for calculation of a far-field')
                error
            end
            
            Coord = zeros(2, N_x*N_y);
            for i=1:N_y
                Coord(1,(i-1)*N_x+1:i*N_x) = x;
            end
            for i=1:N_y
                Coord(2,(i-1)*N_x+1:i*N_x) = y(i);
            end

            obj.v_n = zeros(N_x, N_y);

            N = round(2*obj.R/dr);
            v_n0 = zeros(N,N);
            %phase = phase';

            if strcmp(type,"Dynamic") == 1 % Calculation of a radiation pattern

                for i=1:N
                    for j=1:N
                        if dr*sqrt((i-N/2)^2+(j-N/2)^2)<obj.R
                            v_n0(i,j) = 2*pi*obj.Amplitude_Response.dynamic.f*interp1(obj.r(ceil(obj.N_r/2):end), obj.Amplitude_Response.dynamic.amp(ceil(obj.N_r/2):end,1,2), dr*sqrt((i-N/2)^2+(j-N/2)^2) );
                        end
                    end
                end

                for i=1:Nx
                    for j=1:Ny
                        [~,center_x] = min(abs(x - (x(1)+(i-1)*d+obj.R)));
                        [~,center_y] = min(abs(y - (y(1)+(j-1)*d+obj.R)));
                        obj.v_n(center_x-N/2:center_x+N/2-1, center_y-N/2:center_y+N/2-1) = v_n0*exp(1i*phase(i,j));
                    end
                end

                center = [0.458, 0.519];
                limit1 = 0.2;
                limit2 = 0.8;
                
                normal = direction/sqrt(direction(1)^2+direction(2)^2); % Calculation of a normal vector to the measurement surface

                normal1 = [normal(2), normal(1)];
                arrowStart1 = center;
                arrowEnd1 = center + normal1/10;
                arrow = [arrowStart1; arrowEnd1]';
                arrowStart1 = arrow(1,:);
                arrowEnd1 = arrow(2,:);
                
                n2 = sqrt(1/((normal(2)/normal(1))^2+1));
                n1 = -normal(2)*n2/normal(1);
                if isnan(n1)
                    n1 = 1;
                end
                normal2 = [n2, n1];
                if direction(1)>=0
                    normal2 = -normal2;
                end
                
                arrowStart2 = center + normal2/4;
                for i=1:2
                    if arrowStart2(i)>limit2
                        arrowStart2(i) = limit2;
                    elseif arrowStart2(i) < limit1
                        arrowStart2(i) = limit1;
                    end
                end
                
                arrowEnd2 = center - normal2/4;
                for i=1:2
                    if arrowEnd2(i)>limit2
                        arrowEnd2(i) = limit2;
                    elseif arrowEnd2(i) < limit1
                        arrowEnd2(i) = limit1;
                    end
                end
                arrow = [arrowStart2; arrowEnd2]';
                arrowStart2 = arrow(1,:);
                arrowEnd2 = arrow(2,:);

                figure
                set(gcf, 'Position', [100, 100, 500, 500]);
                imagesc(1e6*y, 1e6*x, real(obj.v_n))
                colorbar
                title('Array of cells (normal velocity, m/s)')
                xlabel('y-coordinate, us')
                ylabel('x-coordinate, us')
                axis image
                shading interp
                annotation('textarrow', arrowStart1, arrowEnd1, 'LineWidth', 2, 'String','normal vector');
                annotation('textarrow', arrowStart2, arrowEnd2, 'Headstyle', 'hypocycloid', 'LineWidth', 1, 'Color', 'r', 'String','Measurement plane');
                f = gca;
                %exportgraphics(f,'Array.png','Resolution',300)
                    
                lambda = c/obj.Amplitude_Response.dynamic.f; % Wavelength
                k = 2*pi/lambda; % Wavenumber

                obj.Frequency_Response.RP.phi = linspace(0,180,obj.N_phi);
                obj.Frequency_Response.RP.p = zeros(1, obj.N_phi);

                if direction(2)<0
                    psi = -acos(normal(1))/pi*180;
                else
                    psi = acos(normal(1))/pi*180;
                end
            
                Plane_coord = zeros(3, obj.N_phi); % Calculation of coordinates on the measurement surface
                for i=1:obj.N_phi
                    Plane_coord(1,i) = -r0*cosd(obj.Frequency_Response.RP.phi(i))*sind(psi);
                    Plane_coord(2,i) = r0*cosd(obj.Frequency_Response.RP.phi(i))*cosd(psi);
                    Plane_coord(3,i) = r0*sind(obj.Frequency_Response.RP.phi(i));
                end

                for n=1:obj.N_phi
                    distance = sqrt( (Coord(1,:)-Plane_coord(1,n)).^2 + (Coord(2,:)-Plane_coord(2,n)).^2 + Plane_coord(3,n).^2 )';
                    obj.Frequency_Response.RP.p(n) =  -1/(2*pi)*sum(obj.v_n(:).*exp(1i*k*distance)./distance)*dr^2 * (-obj.rho_medium*2*pi*obj.Amplitude_Response.dynamic.f);
                end


            elseif strcmp(type,"Sweep") == 1 % Calculation of a frequency response

                obj.Frequency_Response.r = r0;
                obj.Frequency_Response.f = obj.Amplitude_Response.sweep.f;
                obj.Frequency_Response.FR = zeros(1,length(obj.Amplitude_Response.sweep.f));

                for l = 1:length(obj.Amplitude_Response.sweep.f)
                    for i=1:N
                        for j=1:N
                            if dr*sqrt((i-N/2)^2+(j-N/2)^2)<obj.R
                                v_n0(i,j) = 2*pi*obj.Amplitude_Response.sweep.f(l)*interp1(obj.r(ceil(obj.N_r/2):end), obj.Amplitude_Response.sweep.amp(ceil(obj.N_r/2):end,1,l,2), dr*sqrt((i-N/2)^2+(j-N/2)^2) );
                            end
                        end
                    end
    
                    for i=1:Nx
                        for j=1:Ny
                            [~,center_x] = min(abs(x - (x(1)+(i-1)*d+obj.R)));
                            [~,center_y] = min(abs(y - (y(1)+(j-1)*d+obj.R)));
                            obj.v_n(center_x-N/2:center_x+N/2-1, center_y-N/2:center_y+N/2-1) = v_n0;
                        end
                    end
             
                    lambda = c/obj.Amplitude_Response.sweep.f(l); % Wavelength
                    k = 2*pi/lambda; % Wavenumber
    
                    p = -(exp(1i*k*r0)./(2*pi*r0))*trapz(y, trapz(x, obj.v_n, 1)) * (-obj.rho_medium*2*pi*obj.Amplitude_Response.sweep.f(l)); % Pressure at a specified distance

                    obj.Frequency_Response.FR(l) = 10*log10((0.5*p*conj(p))/((p_ref)^2));
                end                
            end
        end


        function obj = Plot_FR(obj, type)

            if strcmp(type,"Dynamic") == 1
            
                if obj.rho_medium == 1.2754
                    p_ref = 20e-6;
                elseif obj.rho_medium == 997
                    p_ref = 1e-6;
                end

                figure
                polarplot(obj.Frequency_Response.RP.phi/180*pi,  10*log10((0.5*obj.Frequency_Response.RP.p.*conj(obj.Frequency_Response.RP.p))/p_ref^2));
                thetalim([0; 180])
                title(['Radiation pattern, ' num2str(1e-6*obj.Amplitude_Response.dynamic.f),' MHz'])
                disp(['Maximum pressure level = ', num2str(10*log10( (0.5*max(obj.Frequency_Response.RP.p)).*conj(max(obj.Frequency_Response.RP.p)) /p_ref^2)), ' dB'])
                f = gca;
                %exportgraphics(f,'DynamicFR.png','Resolution',300)


            elseif strcmp(type,"Sweep") == 1
                
                figure
                plot(1e-6*obj.Frequency_Response.f, real(obj.Frequency_Response.FR))
                title('Frequency response')
                xlabel('Frequency, MHz')
                ylabel('Pressure level, dB')
                xlim(1e-6*[obj.Frequency_Response.f(1); obj.Frequency_Response.f(end)])
                f = gca;
                %exportgraphics(f,'SweepFR.png','Resolution',300)
            end
        end
    end
end