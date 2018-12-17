classdef PointCloudFitter
    properties
        D;  % d = D * u + d0.
        S;  % Shape function matrix.
        d0;
        Hu; % Hu = D' * D; of size nu x nu, with nu = length(nu);
        % We minimize:
        % 0.5 * ||d(u) - y||^2 = 0.5 * u'*Hu*u + fu'*u + ||d0-y||^2
        H;
        K; % Stiffness matrix.
        node_areas; % Array of areas.
        Gp; % d(pv)/du, with pv the internal pressure.
        Aeq;  % Euality contraints.
        Ane;  % Inequality constraints.
        
        % Mesh:
        p_BP0;
        tris;
        normalP0_B;
        
        p_BC; % Camera C position in the bubble frame B.
        
        rhat_B; % Camera ray directions.
        
        % Problem sizes.
        npoints, nbcs, nrays;
        
        % Point cloud standard deviationa and scales.
        sigma_dist;
        T0;
        a;
        cost_factor;
    end   
    
    methods
                   
        % Constructor.
        function this = PointCloudFitter(varargin)
            if (nargin == 16)
                this = this.MakePointCloudFitterFromSeparateData(...
                    varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5},...
                    varargin{6}, varargin{7}, varargin{8}, varargin{9}, varargin{10},...
                    varargin{11}, varargin{12}, varargin{13}, varargin{14}, varargin{15}, varargin{16});
            elseif(nargin == 1)
                file_name = varargin{1};
                fitter_data = load(file_name);
                this = fitter_data.fitter;
            else
                error('Wrong number of arguments.');
            end
        end
        
        % p_BC: position of the camera frame C in the bubble frame B.
        function this = MakePointCloudFitterFromSeparateData(this, ...
                p_BP0, normalP_B, tris, node_boundary, ...
                rhat_C, d0, bar_coos, ray_tri_index, ...
                K, node_areas, Gv, Gp, ...
                sigma_dist, T0, a, p_BC) % scales
            npoints = size(p_BP0, 1);
            nrays = size(rhat_C, 1);
            
            this.p_BP0 = p_BP0;
            this.p_BC = p_BC;
            this.normalP0_B = normalP_B;
            this.sigma_dist = sigma_dist;
            this.T0 = T0;
            this.a = a;
            this.rhat_B = rhat_C;
            this.tris = tris;
            
            % Triplets of values to build D.
            ii = zeros(3*nrays, 1);
            jj = zeros(3*nrays, 1);
            vv = zeros(3*nrays, 1);
            Sv = zeros(3*nrays, 1);  % Shape function matrix.
            
            for iray = 1:nrays
                % Triangle index of the triangle hit by the ray.
                itri = ray_tri_index(iray);
    
                % Node indexes of the triangle hit by the ray.
                tri = tris(itri, :);
    
                % Barycentric coordinates of the intersection point.
                coos2 = bar_coos(:, iray);
                Siray = [1.0 - coos2(1) - coos2(2); coos2(1); coos2(2)];
    
                % Normal interpolated at the iray-th ray.
                normal_iray_B = Siray' * normalP_B(tri, :);
                
                normal_dot_ray = dot(normal_iray_B, rhat_C(iray, :));
                
                %D(iray, tri) = normal_dot_ray * Siray;
                ivalue = 3 * (iray-1) + (1:3);
                ii(ivalue) = [iray, iray, iray];
                jj(ivalue) = tri;
                vv(ivalue) = normal_dot_ray * Siray;
                Sv(ivalue) = Siray;
            end
            
            % This is the factor when u is dimensionless.
            this.cost_factor = 1/nrays * (a/sigma_dist)^2;
            
            this.D = sparse(ii, jj, vv, nrays, npoints, 3*nrays);
            this.S = sparse(ii, jj, Sv, nrays, npoints, 3*nrays);
            this.d0 = d0/a;
            %normalize with the standard deviation and number of rays.
            this.Hu = this.D' * this.D * this.cost_factor;
            
            %       Hu block + pc identity block + pv scalar block.
            Hnnz = nnz(this.Hu) + npoints + 1;
            this.H = spalloc(2*npoints+1, 2*npoints+1, Hnnz);
            this.H(1:npoints, 1:npoints) = this.Hu;
            
            % Regularization.
            %this.H(npoints+(1:npoints), npoints+(1:npoints)) = sparse(1:npoints, 1:npoints, ones(npoints, 1), npoints, npoints);
            %this.H(2*npoints+1, 2*npoints+1) = 0.1;  % minimum pv
            
            
            K = K/T0;
            node_areas = node_areas/(a*a);
            Gp = Gp*a*a/T0;
            
            this.K = K;
            this.node_areas = node_areas;
            this.Gp = Gp;
            
            % Boundary conditions
            ibcs = find(node_boundary > 0);
            nbcs = length(ibcs);
            Abc = sparse(nbcs, npoints);
            for ibc = 1:nbcs
                Abc(ibc, ibcs(ibc)) = 1.0;
            end
            
            % Equality constraints, the actual bubble model.
            % Momentum balance.
            Aeq = sparse(npoints + nbcs + 1, 2*npoints + 1);
            Aeq(1:npoints, 1:npoints) = K;
            %Modify PDEs at BCs            
            %Aeq(1:npoints, ibcs) = 0; % Move known (pc=0) pressure to the RHS.
            %Aeq(ibcs, 1:npoints) = 0; % Replace PDE by BC.
            %for ib=1:nbcs
            %    ig = ibcs(ib);
            %    Aeq(ig, ig) = 1.0;
            %end
            
            % Modify eqs for pc so that it doesnt even show up at BC
            % equation.
            areas_pc = node_areas;
            areas_pc(ibcs) = 0;
            Aeq(1:npoints, npoints+(1:npoints)) = sparse(1:npoints, 1:npoints, -areas_pc, npoints, npoints);
            Aeq(1:npoints, (2*npoints+1)) = -areas_pc;
            
            % BCs on pc, otherwise pc at BCs would just float.
            Aeq(npoints+(1:nbcs), npoints+(1:npoints)) = Abc;
            
            % Volumetric pressure change.
            Gp(ibcs) = 0;
            Gv(ibcs) = 0;
            Aeq(npoints+ nbcs + 1,1:npoints) = Gp';
            Aeq(npoints+ nbcs + 1,2*npoints+1) = -1;
            
            this.Aeq = Aeq;
            
            % Inequality constraint matrix.
            this.Ane = sparse(npoints, 2*npoints + 1);
            this.Ane(1:npoints, (npoints+1):(2*npoints)) = sparse(1:npoints, 1:npoints, ones(1,npoints), npoints, npoints);
            %this.Ane(1:npoints, npoints+ibcs) = 0;
            
            this.npoints = npoints;
            this.nbcs = nbcs;
            this.nrays = nrays;
        end

        
        % y: measured distance.
        function [u, pc, pv, p_BP] = FitPointCloud(this, y)
            npoints = this.npoints;
            nbcs = this.nbcs;
            nrays = this.nrays;
            sigma_dist = this.sigma_dist;            
            a = this.a;
            T0 = this.T0;
                                    
            fu = this.D' * (this.d0 - y/a); % I made fu dimensionless.
            fu = fu * this.cost_factor;
            f = zeros(2*npoints+1, 1);
            f(1:npoints) = fu;
            
            bne = zeros(npoints, 1);
            beq = zeros(npoints + nbcs + 1, 1);
            
            % Make everything dimensionless.
            
            [x, fval, exit_flag] = quadprog(this.H, f, this.Ane, bne, this.Aeq, beq);           
            %[x, fval, exit_flag] = quadprog(this.H, f, [], [], this.Aeq, beq);           
            
            u = x(1:npoints) * a;
            pc = x((npoints+1):(2*npoints)) * T0/a;
            pv = x(2*npoints+1) * T0/a;
            
            p_BP = zeros(npoints, 3);
            for i=1:npoints
                p_BP(i, :) = this.p_BP0(i, :) + u(i) * this.normalP0_B(i, :); 
            end
        end
        
        function [f, fq, fl] = CalcCostFunction(this, y, u, pc, pv)
            npoints = this.npoints;
            nbcs = this.nbcs;
            nrays = this.nrays;
            sigma_dist = this.sigma_dist;            
            
            u = u / this.a;
            y = y / this.a;
            
            fu = this.D' * (this.d0 - y);
            
            fq = 0.5 * u' * this.Hu * u * this.cost_factor; 
            
            fl = fu' * u * this.cost_factor;
            
            f = fq + fl;
        end

        function ui = InterpolateOnPointCloud(this, u)
            % u: quantity to be interpolated. Of size npoints.
            % ui: interplated quantity. Of size nrays.
            ui = this.S * u;
        end
        
    
    end % End of methods section.
    
end