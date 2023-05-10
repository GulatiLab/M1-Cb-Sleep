function [snapshots_n, snapshots_c] = normalize(snapshots)
% Normalize LFP snapshots by subtracting the mean and standard deviation  
% (across trials). Also create a variant where the median is also subracted
% after normalizing.
%
% Inputs:
%       snapshots       - LFP snapshots to be normalized
%
% Outputs:
%       snapshots_n     - Normaized snapshots
%       snapshots_c     - Normaized snapshots with the median removed

    snapshots_n=snapshots;
    
    for j=1:size(snapshots_n,1)
        tmp=reshape(squeeze(snapshots_n(j,:,:)),1,[]);
        snapshots_n(j,:,:)=(snapshots_n(j,:,:)-mean(tmp))./std(tmp);
    end
    
    med=median(snapshots_n,1);
    snapshots_c=bsxfun(@minus,snapshots_n,med);
end