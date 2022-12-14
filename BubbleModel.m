classdef BubbleModel
    properties
        p_BP0;  % Undeformed bubble positions in bubble frame.
        tris;   % Triangle connectivities.
        T0;     % Bubble tension in N/m.
        normal0_B, K, Aeq, node_areas, node_boundary;
        dVdu; % Derivative of the enclosed volume with respect to u.
        V0, p0;
        areas_wbcs; % Node areas vector. At BC nodes they are made equal to zero.
    end
       
    methods
    
        % Constructor.
        % p_BC: position of the camera frame C in the bubble frame B.
        function this = BubbleModel(p_BP0, tris, T0, V0, p0)            
           this.p_BP0 = p_BP0;
           this.tris = tris;
           this.T0 = T0;
           
           if nargin== 5
            this.V0 = V0;
            this.p0 = p0;
           else
            % Thes values are not used for problems wihtout contact.
            % We set them to NaN so that if used by mistake, they lead to a
            % trail of NaNs we can track.
            this.V0 = nan;
            this.p0 = nan;
           end
           
           % Pre-compute quantities that do not change.
           [this.normal0_B, this.K, this.Aeq, this.node_areas, this.node_boundary, this.dVdu, this.areas_wbcs] = ...
               QP_membrane_solver_preproc(p_BP0, tris, T0, this.V0, this.p0);                     
        end                
        
        % Computes deformed bubble solution for a given "fixed" value of
        % internal pressure pv.
        % Solves Eq. 5 in the paper (with pc = 0, no contact).
        function [u, p_BP] = ComputeSolutionWithoutContact(this, pv)
            % Computes displacements, solves Eq. 5.
            rhs = pv * this.areas_wbcs;
            u = this.K\rhs;
            
            % Solve for deformed positions.
            nnodes = length(u);
            p_BP = this.p_BP0;
            for inode = 1:nnodes
                p_BP(inode, :) = p_BP(inode, :) + u(inode) * this.normal0_B(inode, :);
            end
        end
        
        % phi0: Distances from bubble to object. Either negative or close
        % to being negative (object might contact after deformation).
        % H: -dphi/du, a sparse matrix with as many rows as entries in phi0
        % and as many columns as number of unkknows (number of mesh nodes).
        function [u, pr, pv, p_BP, Hmean, lambda] = ComputeSolution(this, phi0, H, Hj)
            
            pa = 101325; % Atmospheric pressure in Pa
            
            %pv = 1.075582560000000e+03; % pv = 0.53 psi - p0
            max_iters = 50;
            relative_tolerance = 1.0e-4;
            omega = 0.7;  % Relaxation. Found by trial and error. orig: 0.7
            pv = 0;
            
            this.V0
            this.p0
            
            % External fixed-point iteration to find pressure.
            for it = 1:max_iters
            
                [u, pv, p_BP, lambda] = QP_membrane_solver(...
                    this.p_BP0, this.tris, this.T0, ...
                    this.normal0_B, this.K, this.Aeq, ...
                    phi0, H, pv, this.areas_wbcs);
                
                % Volume
                V = this.dVdu' * u + this.V0;
                
                % Compute new pressure. Gas ideal law.
                % Use absolute pressure, i.e. add atmospheric pressure.
                P0 = (this.p0+pa);
                pvk = P0 * this.V0 / V - P0;
                
                % Limit pvk, specially for first few iterations.
                pvk = min(pvk, 3*P0);
                pvk = max(0, pvk);  % Pressure should in crease.
                
                %it
                %V
                %pvk
                
                % Convergence error.
                pv_err = pvk - pv;
                pv_rel_err = abs(pv_err/pvk);
                
                %pv_rel_err
                
                % Update pressure. Maybe use relaxation.
                pv = omega * pvk + (1-omega) * pv;
                break;
                if (pv_rel_err < relative_tolerance)
                    sprintf('Converged at iteration %d. pv = %g', it, pv)
                    break;
                end            
            
            end
            
%             if (pv_rel_err > relative_tolerance)
%                 error('DIVERGENCE!!! max-iters: %d. pv = %g', it, pv);
%             end            
            
            npoints = size(this.p_BP0, 1);
            %u = x(1:npoints);
            %pv = x(npoints+1);
            pi = lambda.ineqlin; % KKT multipliers (sparse).
            pr = zeros(npoints, 1);
            pr(Hj) = pi./ this.node_areas(Hj);
            
            %F = this.K * u; % rhs.
            %pr = F ./ this.node_areas; % Actually F = M * p, with M the mass matrix.
            
            % Estimate mean curvature from Lapace's equation: Δp = 2⋅T0⋅H.
            % WARNING: This is not quite the curvature, since when u = 0 then H = 0
            % which is not true.
            Hmean = pr/this.T0/2;
        end
        
        
    
    end % End of methods section.
    
end