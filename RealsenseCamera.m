classdef RealsenseCamera
    properties
        rhat_C;  % Camera rays's directions. nrays x 3.
        nr;  % Total number of rays (pixels). = size(1, rhat_C).
        p_BC_list; % = repmat(p_BC', 1, nr);  % As needed by OPCODE
        
        p_BC;
        tree;
    end
    
    methods
    
        % Constructor.
        function this = RealsenseCamera(varargin)
            if nargin == 3
                this = this.RealsenseCameraFromDataSheet(varargin{1}, varargin{2}, varargin{3});
            elseif nargin == 4
                this = this.RealsenseCameraFromProvidedDirections(varargin{1}, varargin{2}, varargin{3}, varargin{4});
            else
                error('Wrong number of arguments.');
            end
        end            
        
        % p_BC: position of the camera frame C in the bubble frame B.
        function this = RealsenseCameraFromDataSheet(this, p_BC, p_BP, t)
            % Position of the picoflex camera frame C in the bubble frame B.
            this.p_BC = p_BC;
            
            this.rhat_C = generate_realsense_rays();
            this.nr = size(this.rhat_C, 1);
            this.p_BC_list = repmat(p_BC', 1, this.nr);  % As needed by OPCODE

            % Generate AABB tree for the mesh.
            addpath('../opcodemesh/matlab'); % Make OPCODE lib available.
            this.tree = opcodemesh(p_BP', t'); % NOTE!: I am using the transpose!
        end
        
        function this = RealsenseCameraFromProvidedDirections(this, p_BC, rhat_C, p_BP, t)
            % Position of the picoflex camera frame C in the bubble frame B.
            this.p_BC = p_BC;
            
            this.rhat_C = rhat_C; % OPCDE requires the transpose.
            this.nr = size(this.rhat_C, 1);
            this.p_BC_list = repmat(p_BC, 1, this.nr);  % As needed by OPCODE

            % Generate AABB tree for the mesh.
            addpath('../opcodemesh/matlab'); % Make OPCODE lib available.
            this.tree = opcodemesh(p_BP', t'); % NOTE!: I am using the transpose!
        end
        
        function [does_hit, dist, tri_index, bar_coos, p_BY] = GeneratePointCloud(this, p_BP, sigma_dist, default_dist)
            
            this.tree.update(p_BP');
            [does_hit, dist, tri_index, bar_coos, p_BY] = ...
                this.tree.intersect(this.p_BC_list, this.rhat_C');
            p_BY = p_BY';
            % Apparently OPCODE uses single precision. We convert to double
            % since when multiplying by sparse matrices Matlab needs the
            % type to be "double".
            dist = double(dist); 
            
            % Sometimes it seems like OPCODE reports misses (nans). DEBUG!
            nan_indexes = find(isnan(dist));
            dist(nan_indexes) = default_dist(nan_indexes); %this line is the problem
            p_BY(nan_indexes, :) = zeros(length(nan_indexes), 3);

            % Add noise.
            %dist_noise = normrnd(0.0, sigma_dist, nr, 1); % You need a toolbox for
            %this! agghhh!
            dd = zeros(this.nr, 1);
            for ir = 1:this.nr
                dd(ir) = generate_gaussian_samples(0.0, sigma_dist);
                p_BY(ir, :) = p_BY(ir, :) + dd(ir) * this.rhat_C(ir, :);
            end
                        
        end
    
    end % End of methods section.
    
end