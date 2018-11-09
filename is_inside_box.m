function is_inside = is_inside_box(lengths, p_BP)
% Determine if point P, expressed in the box frame B, is inside the box.

half_lengths = lengths / 2;

is_inside = true;

for i = 1:3

if (-half_lengths(i) > p_BP(i) || p_BP(i) > half_lengths(i)) 
    is_inside = false;
    return;
end

end

end