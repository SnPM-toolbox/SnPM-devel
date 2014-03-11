function h=snpm_abline(b,a,varargin)
% FORMAT h = snpm_abline(b,a,...)
% Plots y=a*x+b in dotted line on current axis
% 
% b    - intercept (defaults to 0)
% a    - slope (defaults to 1)
% ...  - Other line options
%
% Like Splus' abline
%
% FORMAT h = snpm_abline('h',y,...)
% Plots horizontal dotted line at y
% 
% FORMAT h = snpm_abline('v',x,...)
% Plots verticle dotted line at x
% 
%
% @(#)snpm_abline.m	1.2 T. Nichols 04/07/08
%$Id: snpm_abline.m,v 8.1 2009/01/29 15:02:57 nichols Exp $	

if (nargin==2) & isstr(b)
  b = lower(b);
else

  if (nargin<1)
    b = 0;
  end
  if (nargin<2)
    a = 0;
  end
  
end

XX=get(gca,'Xlim');
YY=get(gca,'Ylim');

if isstr(b) & (b=='h')

  g=line(XX,[a a],'LineStyle',':',varargin{:});

elseif isstr(b) & (b=='v')

  g=line([a a],YY,'LineStyle',':',varargin{:});

else

  g=line(XX,a*XX+b,'LineStyle',':',varargin{:});

end

if (nargout>0)
  h=g;
end
