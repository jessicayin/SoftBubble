function [phi, nabla_phi] = calc_box_level_set(lengths, p_BP_list)

half_lengths = lengths/2;

npoints = size(p_BP_list, 1);
nabla_phi = zeros(npoints, 3);
phi = -9999 * ones(npoints, 1);

for ipoint = 1:npoints
    p_BP = p_BP_list(ipoint, :)';
        
    for idim = 1:3
        phip = p_BP(idim) - half_lengths(idim);
        phim = -p_BP(idim) - half_lengths(idim);              
        
        if (phi(ipoint) < phip || phi(ipoint) < phim)            
            nabla_phi(ipoint, :) = 0;
            if (phip > phim)
                phi(ipoint) = phip;                
                nabla_phi(ipoint, idim) = 1;
            else
                phi(ipoint) = phim;
                nabla_phi(ipoint, idim) = -1;
            end            
        end      
                
    end
    
end
