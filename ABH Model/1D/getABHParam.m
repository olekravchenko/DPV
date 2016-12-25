function [td, ti, cd, a, b, alpha0, mu] = getABHParam(id)

switch id
    case '01'           
        td      = 1e-3;                                 % T_dT_i / T_e^2Z_d ratio
        ti      = 0.125;                                % T_i/T_e ratio
        cd      = sqrt(td);                             % Dust sound velocity
        a       = 7.5;                                  % Constant in ion-drag force
        b       = 1.6;                                  % Constant in ion-drag force
        alpha0  = 2.0;                                  % Normalized dust-neutral collision frequency
        mu      = 1.5;                                  % Constant mobility
    case '02'           
        td      = 1e-3;                                 % T_dT_i / T_e^2Z_d ratio
        ti      = 0.125;                                % T_i/T_e ratio
        cd      = sqrt(td);                             % Dust sound velocity
        a       = 7.5;                                  % Constant in ion-drag force
        b       = 2.0;                                  % Constant in ion-drag force
        alpha0  = 2.0;                                  % Normalized dust-neutral collision frequency
        mu      = 1.5;                                  % Constant mobility
    otherwise
        td      = 1e-3;                                 % T_dT_i / T_e^2Z_d ratio
        ti      = 0.125;                                % T_i/T_e ratio
        cd      = sqrt(td);                             % Dust sound velocity
        a       = 7.5;                                  % Constant in ion-drag force
        b       = 1.6;                                  % Constant in ion-drag force
        alpha0  = 2.0;                                  % Normalized dust-neutral collision frequency
        mu      = 1.5;                                  % Constant mobility
end