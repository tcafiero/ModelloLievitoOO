classdef YeastClassSpecial < YeastClass
    properties
        %process parameters
        M_f=4.14; c_f=500; y_r=0.5; c=100;
        %Glycolisis_Paramters
        v_g=5.8; k_g=0.27; eta_g=0.64;
        %Fermentation Parameters
        v_f=6.57; k_f=0.16; eta_fe=0.61; eta_fp=0.10; E_max=100;
        %Respiration Parameters
        v_rp=0.83; k_rp=0.18; eta_rp=0.73; v_re=0.11; k_re=0.15; eta_re=0.80;
        %Ethanol
        sigma_e=1.4;
        %Pyruvate
        P_max=1;
        %Inhibitor
        rho=0.02; sigma_i=1.68; I_max=1;                        
        %Metabolic overflow
        a1=0.0002; b1=30;                       
        %Glucose Effect
        a2=0.0002; b2=30;                        
        %Death
        delta=0.1; tau=0.6;                       
        % Reverves
        v_a=0.16; k_a=0.03; eta_a=0.2; a3=0.0002; b3=30; rr=0.3;     
    end
    methods
%       function obj=YeastClassSpecial()
%           obj=obj@YeastClass();
%           self=obj;
%       end
       
    end
    
end
