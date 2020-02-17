function [ d ] = compare( EEG_in, EEG_out )

	threshold = 0.000001;
	
	[numChan, numSamp] = size(EEG_in);
	
	numDiff = length(find(abs(EEG_in-EEG_out) > threshold));
	if (numDiff == 0)
		d = 0;
		return;
	end
	
	numDiff = length(find(abs(EEG_in(:,2:end)-EEG_out(:,1:end-1)) > threshold));
	if (numDiff == 0)
		d = 0;
		return;
	end
	
	numDiff = length(find(abs(EEG_in(:,1:end-1)-EEG_out(:,2:end)) > threshold));
	if (numDiff == 0)
		d = 0;
		return;
	end
	
	d = 1;

end

