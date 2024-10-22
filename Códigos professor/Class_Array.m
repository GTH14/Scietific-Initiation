classdef Class_Array
    properties
        xt  %Positions of each transducer (x-component) [m]
        yt  %Positions of each transducer (y-component) [m]
        zt  %Positions of each transducer (z-component) [m]
        nx  %unit vector normal to the transducer surface (x-component)
        ny  %unit vector normal to the transducer surface (y-component)
        nz  %unit vector normal to the transducer surface (z-component)
        phases %phases applied to each transducer (from 0 to 2*pi)
        Voltage_amplitude = 7; %Voltage amplitude (Peak value)
        c0 =340  %sound velocity - air [m/s]
        rho0=1.2 %fluid density [kg/m^3]
        cp=900; %sound velocity - particle [m/s]
        rhop=30; % particle density [kg/m^3]
        R = .25e-3; %particle radius [m]
        
        freq=40e3 %frequency [Hz]
        a=4.9e-3 %transducer radius [m]
        %Transd_sens =.732; % transducer sensitivity (Acoustic pressure at 1 m from the transducer surface) [Pa/V] (Obs. voltage amplitude Vp) 
        Transd_sens=.462% 1.6 - free-field correction
    end
    properties (Constant)
        Discretization_factor=1000; %delta_x,delta_y,delta_z = wavelength/Discretization factor;
    end
    methods
        % Draw transducers function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Draw_transducers(obj)
            xt=obj.xt;
            yt=obj.yt;
            zt=obj.zt;
            nx=obj.nx;
            ny=obj.ny;
            nz=obj.nz;
            a=obj.a;
            
            teta=0:.1:2*pi;

            for ii=1:size(xt,1)
                %orthonormal basis
                n=[nx(ii) ny(ii) nz(ii)];
                temp1=n(:).'/norm(n);
                temp2=null(temp1).';
                base=[temp1;temp2];
                ez2=base(1,:);
                ey2=base(2,:);
                ex2=-base(3,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %original base vector
                ex1=[1 0 0];
                ey1=[0 1 0];
                ez1=[0 0 1];
                %%%%%%%%%%%%%%

                x=a*cos(teta);
                y=a*sin(teta);
                z=zt(ii)*linspace(0,0,length(teta));

                %Coordinate transformation
                T=[ex2*ex1' ex2*ey1' ex2*ez1';ey2*ex1' ey2*ey1' ey2*ez1';ez2*ex1' ez2*ey1' ez2*ez1'];  
                %%%%%%%%%%%%%%%%%
                X=(T')*[x;y;z];
                x2=X(1,:)+xt(ii);
                y2=X(2,:)+yt(ii);
                z2=X(3,:) + zt(ii)*ones(1,length(teta));
                patch(x2*1000,y2*1000,z2*1000,'r')
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Calculate pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function pressure=Calculate_pressure(obj,xp,yp,zp)
            % xp, yp, zp -> points to evaluate the acoustic pressure
            
            xt=obj.xt;
            yt=obj.yt;
            zt=obj.zt;
            nx=obj.nx;
            ny=obj.ny;
            nz=obj.nz;
            a=obj.a;
            c0 = obj.c0;
            rho0 = obj.rho0;
            Transd_sens = obj.Transd_sens;
            Voltage_amplitude = obj.Voltage_amplitude; 
            phases = obj.phases; 
            
            freq=obj.freq;
            lambda=c0/freq;
            k=2*pi/lambda;
            omega=2*pi*freq;
                        
            N=length(xp); %number of points to evaluate pressure
            for ii=1:N
                r=sqrt((xp(ii)-xt).^2 + (yp(ii)-yt).^2 + (zp(ii)-zt).^2);
                teta=acos( (nx.*(xp(ii)-xt) + ny.*(yp(ii)-yt) + nz.*(zp(ii)-zt) )./r);
                pressure(ii)=sum(Voltage_amplitude.*(Transd_sens./r).*(2*besselj(1,k*a*sin(teta+eps))./(k*a*sin(teta+eps))).*exp(+j*(k*r - phases)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Calculate particle velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ux,uy,uz]=Calculate_velocity(obj,xp,yp,zp)
            % xp, yp, zp -> points to evaluate the particle velocity
            
            % ux, uy, uz -> particle velocity in the x, y and z directions
        
            xt=obj.xt;
            yt=obj.yt;
            zt=obj.zt;
            nx=obj.nx;
            ny=obj.ny;
            nz=obj.nz;
            a=obj.a;
            c0 = obj.c0;
            rho0 = obj.rho0;
            Transd_sens = obj.Transd_sens;
            Voltage_amplitude = obj.Voltage_amplitude; 
            phases = obj.phases; 
            Discretization_factor = obj.Discretization_factor; 
            
            freq=obj.freq;
            lambda=c0/freq;
            k=2*pi/lambda;
            omega=2*pi*freq;
            
            N=length(xp); %number of points to evaluate the particle velocity
            delta_x=lambda/Discretization_factor; 
            delta_y=lambda/Discretization_factor; 
            delta_z=lambda/Discretization_factor; 
            for ii=1:N
                veloc_potential_0=-Calculate_pressure(obj,xp(ii),yp(ii),zp(ii))/(j*omega*rho0);
                veloc_potential_xplus =-Calculate_pressure(obj,xp(ii) + delta_x,yp(ii),zp(ii))/(j*omega*rho0);
                veloc_potential_xminus =-Calculate_pressure(obj,xp(ii) - delta_x,yp(ii),zp(ii))/(j*omega*rho0);
                veloc_potential_yplus =-Calculate_pressure(obj,xp(ii),yp(ii) + delta_y,zp(ii))/(j*omega*rho0);
                veloc_potential_yminus =-Calculate_pressure(obj,xp(ii),yp(ii) - delta_y,zp(ii))/(j*omega*rho0);
                veloc_potential_zplus =-Calculate_pressure(obj,xp(ii),yp(ii),zp(ii) + delta_z)/(j*omega*rho0);
                veloc_potential_zminus =-Calculate_pressure(obj,xp(ii),yp(ii),zp(ii) - delta_z)/(j*omega*rho0);
                
                ux(ii) = (veloc_potential_xplus - veloc_potential_xminus)/(2*delta_x);
                uy(ii) = (veloc_potential_yplus - veloc_potential_yminus)/(2*delta_y);
                uz(ii) = (veloc_potential_zplus - veloc_potential_zminus)/(2*delta_z);
            end
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        % Calculate Gorkov Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function U=Calculate_Gorkov_potential(obj,xp,yp,zp)
            % xp, yp, zp -> points to evaluate the Gorkov potential
            
            xt=obj.xt;
            yt=obj.yt;
            zt=obj.zt;
            nx=obj.nx;
            ny=obj.ny;
            nz=obj.nz;
            a=obj.a;
            c0 = obj.c0;
            rho0 = obj.rho0;
            cp = obj.cp;
            rhop = obj.rhop;
            R = obj.R;
            Transd_sens = obj.Transd_sens;
            Voltage_amplitude = obj.Voltage_amplitude; 
            phases = obj.phases; 
            Discretization_factor = obj.Discretization_factor; 
            
            freq=obj.freq;
            lambda=c0/freq;
            k=2*pi/lambda;
            omega=2*pi*freq;
            
            
            N=length(xp);
            delta_x=lambda/Discretization_factor; 
            delta_y=lambda/Discretization_factor; 
            delta_z=lambda/Discretization_factor; 
            for ii=1:N
                
                pressure=Calculate_pressure(obj,xp(ii),yp(ii),zp(ii));
                [ux,uy,uz]=Calculate_velocity(obj,xp(ii),yp(ii),zp(ii));
            
                p_quad=(abs(pressure).^2)/2;
                %mod_u=sqrt(ux.^2 +uy.^2 + uz.^2);
                %u_quad=(abs(mod_u).^2)/2;
                u_quad=(ux*conj(ux) + uy*conj(uy) + uz*conj(uz))/2;
                
                const=2*pi*(R^3);
                %f1=1; %TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO
                %f2=1; %TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO
                f1 = 1 - (rho0*(c0^2))/(rhop*(cp^2));
                f2 = 2*(rhop-rho0)/(2*rhop+rho0);
                U(ii)=const*( f1*p_quad/(3*rho0*c0^2) - f2*rho0*u_quad/2);
                
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        
        % Calculate Radiation Force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Fx,Fy,Fz]=Calculate_radiation_force_Gorkov(obj,xp,yp,zp)
            % xp, yp, zp -> points to evaluate the particle velocity
            
            % Fx, Fy, Fz -> Fadiation force in the x, y and z directions
            
            xt=obj.xt;
            yt=obj.yt;
            zt=obj.zt;
            nx=obj.nx;
            ny=obj.ny;
            nz=obj.nz;
            a=obj.a;
            c0 = obj.c0;
            rho0 = obj.rho0;
            Transd_sens = obj.Transd_sens;
            Voltage_amplitude = obj.Voltage_amplitude; 
            phases = obj.phases; 
            Discretization_factor = obj.Discretization_factor; 
            
            freq=obj.freq;
            lambda=c0/freq;
            k=2*pi/lambda;
            omega=2*pi*freq;
            
            N=length(xp); %number of points to evaluate the particle velocity
            delta_x=lambda/Discretization_factor; 
            delta_y=lambda/Discretization_factor; 
            delta_z=lambda/Discretization_factor; 
            for ii=1:N
                Fx(ii) = - (Calculate_Gorkov_potential(obj,xp(ii) + delta_x,yp(ii),zp(ii)) - Calculate_Gorkov_potential(obj,xp(ii) - delta_x,yp(ii),zp(ii))     )/(2*delta_x);
                Fy(ii) = - (Calculate_Gorkov_potential(obj,xp(ii),yp(ii)+delta_y,zp(ii)) - Calculate_Gorkov_potential(obj,xp(ii),yp(ii)-delta_y,zp(ii)))/(2*delta_y);
                Fz(ii) = - (Calculate_Gorkov_potential(obj,xp(ii),yp(ii),zp(ii)+delta_z) - Calculate_Gorkov_potential(obj,xp(ii),yp(ii),zp(ii)-delta_z))/(2*delta_z);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Calculate Radiation Force - Acoustokinetics %%%%%%%%%%%%%%%%%%%%%
        function [Fx,Fy,Fz]=Calculate_radiation_force_Acoustokinetics(obj,xp,yp,zp)
            % xp, yp, zp -> points to evaluate the particle velocity
            
            % Fx, Fy, Fz -> Fadiation force in the x, y and z directions
            
            xt=obj.xt;
            yt=obj.yt;
            zt=obj.zt;
            nx=obj.nx;
            ny=obj.ny;
            nz=obj.nz;
            a=obj.a;
            c0 = obj.c0;
            rho0 = obj.rho0;
            cp = obj.cp;
            rhop = obj.rhop;
            R = obj.R;
            Transd_sens = obj.Transd_sens;
            Voltage_amplitude = obj.Voltage_amplitude; 
            phases = obj.phases; 
            Discretization_factor = obj.Discretization_factor; 
            
            %f0=1; %TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO (comparar eq. de Gorkov)
            %f1=1; %TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO TEMPORARIO (comparar eq. de Gorkov)
            f0 = 1 - (rho0*(c0^2))/(rhop*(cp^2));
            f1 = 2*(rhop-rho0)/(2*rhop+rho0);
            
            freq=obj.freq;
            lambda=c0/freq;
            k=2*pi/lambda;
            omega=2*pi*freq;
            
            alpha_a=((4*pi*(R^3))/(3*rho0*(c0^2)))*f0*(-1 + j*(1/3)*(f0+f1)*((k*R)^3));
            beta_a=((2*pi*(R^3))/(rho0*(c0^2)))*f1*(1 + j*(1/6)*(f1)*((k*R)^3)); 
            
            
            N=length(xp); %number of points to evaluate the particle velocity
            delta_x=lambda/Discretization_factor; 
            delta_y=lambda/Discretization_factor; 
            delta_z=lambda/Discretization_factor; 
            for ii=1:N
                pressure_0 = Calculate_pressure(obj,xp(ii),yp(ii),zp(ii));
                pressure_xplus = Calculate_pressure(obj,xp(ii) + delta_x,yp(ii),zp(ii));
                pressure_xminus = Calculate_pressure(obj,xp(ii) - delta_x,yp(ii),zp(ii));
                pressure_yplus = Calculate_pressure(obj,xp(ii),yp(ii) + delta_y,zp(ii));
                pressure_yminus = Calculate_pressure(obj,xp(ii),yp(ii) - delta_y,zp(ii));
                pressure_zplus = Calculate_pressure(obj,xp(ii),yp(ii),zp(ii) + delta_z);
                pressure_zminus = Calculate_pressure(obj,xp(ii),yp(ii),zp(ii) - delta_z);
                
                pressure_xplus_yplus = Calculate_pressure(obj,xp(ii) + delta_x,yp(ii) + delta_y,zp(ii));
                pressure_xplus_zplus = Calculate_pressure(obj,xp(ii) + delta_x,yp(ii) ,zp(ii) + delta_z);
                pressure_xminus_yminus = Calculate_pressure(obj,xp(ii) - delta_x,yp(ii) - delta_y,zp(ii));
                pressure_xminus_zminus = Calculate_pressure(obj,xp(ii) - delta_x,yp(ii) ,zp(ii) - delta_z);
                pressure_yplus_zplus = Calculate_pressure(obj,xp(ii),yp(ii) + delta_y,zp(ii) + delta_z);
                pressure_yminus_zminus = Calculate_pressure(obj,xp(ii),yp(ii) - delta_y,zp(ii) - delta_z);
                
                grad_pconj_x = (conj(pressure_xplus) - conj(pressure_xminus))/(2*delta_x);
                grad_pconj_y = (conj(pressure_yplus) - conj(pressure_yminus))/(2*delta_y);
                grad_pconj_z = (conj(pressure_zplus) - conj(pressure_zminus))/(2*delta_z);
                
                d2_pconj_dx2 = (conj(pressure_xplus) - 2*conj(pressure_0) + conj(pressure_xminus))/(delta_x^2);
                d2_pconj_dy2 = (conj(pressure_yplus) - 2*conj(pressure_0) + conj(pressure_yminus))/(delta_y^2);
                d2_pconj_dz2 = (conj(pressure_zplus) - 2*conj(pressure_0) + conj(pressure_zminus))/(delta_z^2);
                
                d2_pconj_dxdy = (conj(pressure_xplus_yplus) - conj(pressure_xplus) - conj(pressure_yplus)...
                    + 2*conj(pressure_0) -conj(pressure_xminus) - conj(pressure_yminus) + ...
                    conj(pressure_xminus_yminus))/(2*delta_x*delta_y);
                d2_pconj_dxdz = (conj(pressure_xplus_zplus) - conj(pressure_xplus) - conj(pressure_zplus)...
                    + 2*conj(pressure_0) -conj(pressure_xminus) - conj(pressure_zminus) + ...
                    conj(pressure_xminus_zminus))/(2*delta_x*delta_z);
                d2_pconj_dydz = (conj(pressure_yplus_zplus) - conj(pressure_yplus) - conj(pressure_zplus)...
                    + 2*conj(pressure_0) -conj(pressure_yminus) - conj(pressure_zminus) +...
                    conj(pressure_yminus_zminus))/(2*delta_y*delta_z);
                                
                grad_p_x = (pressure_xplus - pressure_xminus)/(2*delta_x);
                grad_p_y = (pressure_yplus - pressure_yminus)/(2*delta_y);
                grad_p_z = (pressure_zplus - pressure_zminus)/(2*delta_z);
                
                first_term_x=alpha_a*pressure_0*grad_pconj_x;
                first_term_y=alpha_a*pressure_0*grad_pconj_y;
                first_term_z=alpha_a*pressure_0*grad_pconj_z;
                

                
                second_term_x=beta_a*(1/(k^2))*(grad_p_x*d2_pconj_dx2  + grad_p_y*d2_pconj_dxdy + grad_p_z*d2_pconj_dxdz  );
                second_term_y=beta_a*(1/(k^2))*(grad_p_x*d2_pconj_dxdy + grad_p_y*d2_pconj_dy2  + grad_p_z*d2_pconj_dydz  );
                second_term_z=beta_a*(1/(k^2))*(grad_p_x*d2_pconj_dxdz + grad_p_y*d2_pconj_dydz + grad_p_z*d2_pconj_dz2  );
                
                Fx(ii)=(1/2)*real(first_term_x + second_term_x);
                Fy(ii)=(1/2)*real(first_term_y + second_term_y);
                Fz(ii)=(1/2)*real(first_term_z + second_term_z);
                                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    
    
end