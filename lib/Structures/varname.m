function out = varname(varargin)
out=[];
    for n=1:nargin
  out{n} = inputname(n);
    end
end