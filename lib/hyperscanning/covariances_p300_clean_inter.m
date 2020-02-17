function Xout=covariances_p300_clean_inter(X,nusers,P300_ref_orientation)
if nargin <3
    P300_ref_orientation='commonP1';
end
if nargin < 2
    nusers=1;
end
switch P300_ref_orientation
    case 'commonP1'
        chan_nb=size(X,1)/(nusers+1);
        Xout=X;
        for Block=1:(nusers)
            example{Block}=((Block)*chan_nb+1):(Block+1)*chan_nb;
        end
        C=combnk(example,2);
        for i=1:size(C,1)
            Xout(C{i,1},C{i,2},:)=0;
            Xout(C{i,2},C{i,1},:)=0;
        end
        disp([ 'Inter Statistics have been set to ZERO for ' P300_ref_orientation])
    case {'multiP1','noP1'}
         chan_nb=size(X,1)/(nusers);
        Xout=X;
        for Block=1:(nusers)
            example{Block}=((Block-1)*chan_nb+1):(Block)*chan_nb;
        end
        C=combnk(example,2);
        for i=1:size(C,1)
            Xout(C{i,1},C{i,2},:)=0;
            Xout(C{i,2},C{i,1},:)=0;
        end
        disp([ 'Inter Statistics have been set to ZERO for ' P300_ref_orientation])
        
    otherwise
        fprintf('covariances_p300_clean_inter : \n Error in P300_ref_orientation \n')
        
end

end