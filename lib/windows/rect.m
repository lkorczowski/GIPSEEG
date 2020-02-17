function xout=rect(varargin)
if nargin==3
    a=varargin{1};
    b=varargin{2};
    x=varargin{3};
elseif nargin==2
    a=1;
    b=varargin{1};
    x=varargin{2};
elseif nargin==1
    a=1;
    b=1;
    x=varargin{1};
end
        
xout=abs(x)<b;
xout=xout*a;