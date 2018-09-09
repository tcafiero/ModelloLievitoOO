classdef YeastClass < handle
	properties (Abstract=true)
        M_f, c_f, y_r, c, v_g, k_g, eta_g, k_re, v_f,...
    	eta_rp, eta_fe, E_max, k_rp, v_rp, k_f, v_re,  sigma_e, P_max, rho, sigma_i,...
    	I_max, a1, b1, a2, b2, delta, tau, eta_re, eta_fp, v_a, k_a, eta_a, a3, b3, rr
	end
	properties %(Access=private)
		B              	%g Viable microbial mass
		M        		%g Total microbial mass
		c_G          	%g/l Glucose Concentration
		c_E          	%g/l Ethanol Concentration
		c_P      		%Pyruvate concentration inside cell microbial mass g/l
		c_I            	%Enviromental toxic concentration g/l
		c_R      		%Reserves concentration
		R_max           %Max amount of reserves
        F0
    end
	methods (Static)
		function obj=YeastClass()
        end
        
		function dx_dt = Lieviti_eqs (self, t1, x1)
			% Intermediate variables
            global t_f mu;

            % Calculate concentrations            
			self.B=x1(4)+x1(5);              %g Viable microbial mass
			self.M=x1(4)+x1(5)+x1(8);        %g Total microbial mass
			self.c_G=(x1(2)/x1(1));          %g/l Glucose Concentration
			self.c_E=(x1(3)/x1(1));          %g/l Ethanol Concentration
			self.c_P=((x1(4)/(self.B+x1(8)))*self.c);      %Pyruvate concentration inside cell microbial mass g/l
			self.c_I=x1(6)/x1(1);            %Enviromental toxic concentration g/l
			self.c_R=x1(8)/(self.B+x1(8))*self.c;      %Reserves concentration
			self.R_max=self.rr*self.M;                %Max amount of reserves

			n_e=abs(self.sigma_e)*self.c_E/abs(self.E_max);
			n_i=abs(self.sigma_i)*self.c_I/abs(self.I_max);        % Inhibitor Negative Feedback
			mo=inv(1+abs(self.a1)*exp(abs(self.b1)*(self.c_P)));      % Metabolic Overflow
			ge=inv(1+abs(self.a2)*exp(abs(self.b2)*self.c_P));      % Glucose Effect
			ra=(1-(inv(1+abs(self.a3)*exp(abs(self.b3)*self.c_P))));  % Reserves accumulation switch
			d=(self.c_P)>self.tau;                            % Death switch
			if t1<=t_f
				self.F0=0;
			else
				self.F0=self.M_f*mu/(self.c_f*self.y_r);
			end

			% MODEL FLUXES
			Feeding=self.c_f*self.F0*exp(mu*(t1-t_f));
			Uptake_G=abs(self.v_g)*self.c_G/(abs(self.k_g)+self.c_G)*self.B*(1-self.c_P/self.P_max)*(1-n_e);
			Respiration_P=abs(self.v_rp)*self.c_P/(abs(self.k_rp)+self.c_P)*self.B*(1-n_e)*(1-n_i)*(ge);
			Fermentation=abs(self.v_f)*self.c_P/(abs(self.k_f)+self.c_P)*self.B*(1-n_e)*(1-n_i)*(1-mo);
			Respiration_E=abs(self.v_re)*self.c_E/(abs(self.k_re)+self.c_E)*self.B*(1-n_e)*(1-n_i)*(ge);
			Secretion=self.rho*(min(0.9,abs(self.eta_rp))*Respiration_P+min(0.9,abs(self.eta_re))*Respiration_E+min(0.9,abs(self.eta_fp))*Fermentation);
			Accumulation=abs(self.v_a)*self.c_P/(abs(self.k_a)+self.c_P)*self.B*(1-x1(8)*inv(self.R_max))*ra;
			Death_P=d*self.delta*x1(4);
			Death_M=d*self.delta*x1(5);
			Death_R=d*self.delta*x1(8);

			% MODEL EQUATIONS
			dV=self.F0*exp(mu*(t1-t_f));
			dG=Feeding-Uptake_G;
			dE=min(0.9,abs(self.eta_fe))*Fermentation-Respiration_E;
			dP=min(0.9,abs(self.eta_g))*Uptake_G-Respiration_P-Fermentation-Accumulation-Death_P;
			dCm=min(0.9,abs(self.eta_rp))*Respiration_P+min(0.9,abs(self.eta_re))*Respiration_E+min(0.9,abs(self.eta_fp))*Fermentation-Secretion-Death_M;
			dI=Secretion;
			dD=Death_P+Death_M+Death_R;
			dR=min(0.9,abs(self.eta_a))*Accumulation-Death_R;
			dx_dt = [dV;dG;dE;dP;dCm;dI;dD;dR];
		end
	end
end
