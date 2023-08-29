classdef Constants

    % System variables
    properties (Constant)
        % Universal Gas Constant in JK^-1mol^-1
        R = 8.314;
        
        % Faradays Constant
        F = 96487;
        
        % Bruggeman factor in different battery regions
        brugn = 1.5;
        brugs = 4.0;
        brugp = 1.5;
        
        % System Variables
        
        
        % System temperature in K
        T = 298;
        
        % Active material fraction in -ve & +ve electrode
        e1n = 0.5;
        e1p = 0.49;
        
        % Volume fraction in different battery regions
        e2n = 0.33;
        e2s = 0.54;
        e2p = 0.332;
        
        % Transferennce number of Li ion species dissolved in liquid
        t_plus = 0.363;
        
        % Characteristic radius of electrode particles in m
        rn = 12.5e-6;
        rp = 8.5e-6;
        
        % Thickness of different electrode regions in m
        ln = 120e-6;
        ls = 30e-6;
        lp = 150e-6;
        
        % Total cell thickness
        L = Constants.ln + Constants.ls + Constants.lp;
        
        % Maximum surface solid phase concentration in mol m^-3
        csn_max=26390;
        csp_max=22860;
        
        % Initial electrolyte concentration in mol m^-3
        c20 = 1200;
        
        % Electronic conductivity 
        k1n = 100 * (Constants.e1n ^ Constants.brugn);
        k1p = 3.8 * (Constants.e1p ^ Constants.brugp);
        
        % Specific surface area of active material in negative electrode
        an = 3 * Constants.e1n / Constants.rn;
        
        % Specific surface area of active material in positive electrode
        ap = 3 * Constants.e1p / Constants.rp;
        
        % Miscellaneous
        theta = Constants.R * Constants.T * (1 - Constants.t_plus)/Constants.F;
    end


    % System functions
    methods (Static)

        % Average reaction rates
        function jnavg = jn_avg(I)
            jnavg = I/(Constants.an*Constants.F*Constants.ln);
        end

        function jpavg = jp_avg(I)
            jpavg = I/(Constants.ap*Constants.F*Constants.lp);
        end

        % Ionic conductivity
        function K2 = k2(c2) 
            K2 = 1e-4 * c2 * ( ...
                -10.5 + 0.074 * Constants.T - 6.69e-5 * Constants.T^2 + ...
                6.68e-4 * c2 -1.78e-5 * c2 * Constants.T + ...
                2.8e-8 * c2 * Constants.T^2 + 4.94e-7 * c2.^2 - ...
                8.86e-10 * c2.^2 * Constants.T ...
            ).^2;
        end

        % Ionic conductivity in different regions of battery
        function K2N = k2n(c2n) 
            K2N= Constants.k2(c2n) * (Constants.e2n ^ Constants.brugn);
        end

        function K2S = k2s(c2s) 
            K2S = Constants.k2(c2s) * (Constants.e2s ^ Constants.brugs);
        end

        function K2P = k2p(c2p) 
            K2P = Constants.k2(c2p) * (Constants.e2p ^ Constants.brugp);
        end

        % Electrolyte diffusivity of material 
        function d2 = D2(c2)
            d2 = 1e-4 * 10 ^ (-2.2e-4 * c2 -4.43 * (54/(Constants.T-229-0.05 * c2)));
        end

        % Electrolyte diffusivity in different regions of battery
        function d2n = D2n(c2n) 
            d2n = Constants.D2(c2n) * (Constants.e2n ^ Constants.brugn);
        end

        function d2s = D2s(c2s) 
            d2s = Constants.D2(c2s) * (Constants.e2s ^ Constants.brugs);
        end

        function d2p = D2p(c2p) 
            d2p = Constants.D2(c2p) * (Constants.e2p ^ Constants.brugp);
        end
        
        % State of charge in electrodes
        function socn = SoCn(csn)
            socn = csn/Constants.csn_max;
        end

        function socp = SoCp(csp)
            socp = csp/Constants.csp_max;
        end

        % Open circuit potentials
        function un = Un(SOCn)
            un = 0.16 + 1.32*exp(-3*SOCn) + 10*exp(-2000*SOCn);
        end

        function up = Up(SOCp)
            up = 4.1983 + 0.0565*tanh(-14.554*SOCp + 8.6094) ...
            - 0.0275*(-1.9011 + 1/(0.9984-SOCp).^0.4924) ...
            -0.1571*exp(-0.0474*SOCp.^8) + 0.8102*exp(-40*(SOCp-0.1339));
        end
    end
end
