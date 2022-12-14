function [p_BP, u, bubble] = inflate_membrane(p_BP0, tris, h0)

% The model is linear. That is, deformations will be proportional to pv/T0.
% Therefore we provide unit values of these parameters to obtain the
% "shape" to then scale it to h0.
T0 = 1.0;
pv = 1.0;


% Construct a model without contact constraints.
bubble = BubbleModel(p_BP0, tris, T0);

% Compute deformed configuration due to pressure pv.
[u, p_BP] = bubble.ComputeSolutionWithoutContact(pv);

z=p_BP(:,3);
z_max = max(z);

% Renormalize.
p_BP(:,3) = z/z_max * h0;
