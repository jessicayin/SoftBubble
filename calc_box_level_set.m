function phi = calc_box_level_set(lengths, p_BP)

half_lengths = lengths/2;

phi = zeros(length(p_BP), 1);
for side=1:3
    pp = p_BP(:, side) - half_lengths(side);
    pm = -p_BP(:, side) - half_lengths(side);

    p = max(pp, pm);
    if (side == 1)
        phi = p;
    else        
        phi = max(phi, p);
    end
end
