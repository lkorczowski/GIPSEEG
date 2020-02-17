function d = SYNC(P1,P2,method)

if (nargin<3)||(isempty(method))
    method_dist = 'euclid';
end
if (nargin<4)
    arg_dist = {};
end

switch method_dist
    case 'riemann'
        d = SYNC_Riem(P1,P2);
    case 'kullback'
        d = distance_kullback(C1,C2);
    case 'logeuclid'
        d = distance_logeuclid(C1,C2);
    case 'opttransp'
        d = distance_opttransp(C1,C2);
    case 'ld'
        d = distance_ld(C1,C2);
    otherwise
        d = sqrt(norm(C1-C2,'fro'));
        disp('WARNING : unknown distance setup... use Frobenius norm')
end