function slope = computeInitialSlope(AoA, Va_I)

y = 1;
x = tan(AoA) * y;

Va_B = [x;y];
Va_B = Va_B ./ norm(Va_B);

Va_I = Va_I./norm(Va_I);

% slope = acos(dot(Va_B, -Va_I));
slope =  sign(Va_I(1))*asin(norm(cross([Va_B;0], [Va_I; 0])));

end