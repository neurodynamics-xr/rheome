function tf = sel_pred(V, op, theta)
% SEL_PRED  V op theta, elementwise.
% Author: Diellor Basha, 2026
    switch op
        case '>',  tf = V > theta;   case '>=', tf = V >= theta;
        case '<',  tf = V < theta;   case '<=', tf = V <= theta;
    end
end
% Author: Diellor Basha, 2026
