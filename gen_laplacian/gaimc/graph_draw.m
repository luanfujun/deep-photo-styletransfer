function h = graph_draw(adj, xy, varargin)
% GRAPH_DRAW Draw a picture of a graph when the coordinates are known
%
% graph_draw(A, xy) draws a picture of graph A where node i is placed
% at x = xy(i,1), y = xy(i,2).  In the drawing, shaded nodes have 
% self loops.
%
% Some of the parameters of the drawing are controlled by specifying
% optional parameters in the call graph_draw(A, xy, key, value).  The keys
% and default values are 
%      'linestyle'   -  default '-' 
%      'linewidth'   -  default .5
%      'linecolor'   -  default Black
%      'fontsize'    -  fontsize for labels, default 8 
%      'labels'      -  Cell array containing labels <Default : '1':'N'>
%      'shapes'      -  1 if node is a box, 0 if oval <Default : zeros>
% 
% h = graph_draw(A,xy,...) returns a handle for each object.  h(i,1) is
% the text handle for vertex i, and h(i,2) is the circle handle for
% vertex i.  
%
% Originally written by Erik A. Johnson, Ali Taylan Cemgil, and Leon Peskin
% Modified by David F. Gleich for gaimc package.
%
% See also GPLOT
% 
% Example:
%   load_gaimc_graph('dfs_example');
%   graph_draw(A,xy);


% 2009-02-26 interface modified by David Gleich <dgleich@stanford.edu> 
%            to remove automatic layout
% 2009-05-15: Added example

% 24 Feb 2004  cleaned up, optimized and corrected by Leon Peshkin pesha @ ai.mit.edu 
% Apr-2000  draw_graph   Ali Taylan Cemgil   <cemgil@mbfys.kun.nl> 
% 1995-1997 arrow        Erik A. Johnson     <johnsone@uiuc.edu>

linestyle = '-';       %   --   -. 
linewidth = .5;        %   2 
linecolor = 'Black';   %   Red
fontsize = 8;
N = size(adj,1);
color = ones(N, 3);         % colors of elipses around text 
labels = cellstr(int2str((1:N)'));    %  labels = cellstr(char(zeros(N,1)+double('+')));
node_t = zeros(N,1);                  %  
for i = 1:2:length(varargin)                  % get optional args
    switch varargin{i}
        case 'linestyle', linestyle = varargin{i+1};
        case 'linewidth', linewidth = varargin{i+1};
        case 'linecolor', linecolor = varargin{i+1};
        case 'labels', labels  = varargin{i+1};
        case 'fontsize',  fontsize = varargin{i+1}; 
        case 'shapes', node_t  = varargin{i+1};  node_t = node_t(:);
    end
end

x = xy(:,1);
x = x - min(x);
y = xy(:,2);
y = y - min(y);
% scale the graph so it's between 0 and 1
xrange = max(x);
yrange = max(y);
scalefactor = max(xrange,yrange);
x = x/scalefactor;
y = y/scalefactor;

lp_ndx = find(diag(adj));       %  recover from self-loops = diagonal ones 
color(lp_ndx,:) = repmat([.8 .8 .8],length(lp_ndx),1);     %  makes self-looped nodes blue
adj = adj - diag(diag(adj));    % clean up the diagonal 

axis([-0.1 1.1 -0.1 1.1]);
axis off; 
set(gcf,'Color',[1 1 1]);
set(gca,'XTick',[], 'YTick',[], 'box','on'); % axis('square');   %colormap(flipud(gray));

idx1 = find(node_t == 0); wd1 = [];   %  Draw  nodes 
if ~isempty(idx1),
    [h1 wd1] = textoval(x(idx1), y(idx1), labels(idx1), fontsize, color);
end;

idx2 = find(node_t ~= 0); wd2 = [];
if ~isempty(idx2),
    [h2 wd2] = textbox(x(idx2), y(idx2), labels(idx2), color);
end;

wd = zeros(size(wd1,1) + size(wd2,1),2);
if ~isempty(idx1), wd(idx1, :) = wd1; end;
if ~isempty(idx2), wd(idx2, :) = wd2; end;

for node = 1:N                          %  Draw  edges 
  edges = find(adj(node,:) == 1);
  for node2 = edges
    sign = 1;
    if ((x(node2) - x(node)) == 0)
	    if (y(node) > y(node2)), alpha = -pi/2; else alpha = pi/2; end;
    else
        alpha = atan((y(node2)-y(node))/(x(node2)-x(node)));
	    if (x(node2) <= x(node)), sign = -1; end;
    end;
    dy1 = sign.*wd(node,2).*sin(alpha);   dx1 = sign.*wd(node,1).*cos(alpha);
    dy2 = sign.*wd(node2,2).*sin(alpha);  dx2 = sign.*wd(node2,1).*cos(alpha);    
    if  (adj(node2,node) == 0)           % if directed edge
        my_arrow([x(node)+dx1 y(node)+dy1], [x(node2)-dx2 y(node2)-dy2]);
    else	   
        line([x(node)+dx1 x(node2)-dx2], [y(node)+dy1 y(node2)-dy2], ...
            'Color', linecolor, 'LineStyle', linestyle, 'LineWidth', linewidth);
        adj(node2,node) = -1;         % Prevent drawing lines twice
    end;
  end;
end;

if nargout > 2
    h = zeros(length(wd),2);
    if ~isempty(idx1), h(idx1,:) = h1;   end;
    if ~isempty(idx2), h(idx2,:) = h2;   end;
end;

function [t, wd] = textoval(x, y, str, fontsize, c)
%  [t, wd] = textoval(x, y, str, fontsize)    Draws an oval around text objects
% INPUT:   x, y - Coordinates
%           str - Strings
%             c - colors 
% OUTPUT:     t - Object Handles
%         width - x and y  width of ovals 
if ~isa(str,'cell'), str = cellstr(str); end;
N = length(str);    
wd = zeros(N,2);
temp = zeros(N,2);
for i = 1:N,
    tx = text(x(i),y(i),str{i},'HorizontalAlignment','center','VerticalAlign','middle','FontSize', fontsize);
    sz = get(tx, 'Extent');
    wy = sz(4);
    wx = max(2/3*sz(3), wy); 
    wx = 0.9 * wx;        %  might want to play with this .9 and .5 coefficients 
    wy = 0.5 * wy;
    ptc = ellipse(x(i), y(i), wx, wy, c(i,:));
    set(ptc, 'FaceColor', c(i,:));   % 'w'
    wd(i,:) = [wx wy];
    delete(tx);
    tx = text(x(i),y(i),str{i},'HorizontalAlignment','center','VerticalAlign','middle', 'FontSize', fontsize);      
    temp(i,:) = [tx ptc];
end;
t = temp; 

function [p] = ellipse(x, y, rx, ry, c)
%  [p] = ellipse(x, y, rx, ry)    Draws Ellipse shaped patch objects
% INPUT:  x,y -  N x 1 vectors of x and y coordinates
%      Rx, Ry -   Radii
%           C -  colors 
% OUTPUT:   p -   Handles of Ellipse shaped path objects

  if length(rx)== 1, rx = ones(size(x)).*rx; end;
  if length(ry)== 1, ry = ones(size(x)).*ry; end;
N = length(x);
p = zeros(size(x));
t = 0:pi/30:2*pi;
for i = 1:N
	px = rx(i) * cos(t) + x(i);    py = ry(i) * sin(t) + y(i);
	p(i) = patch(px, py, c(i,:));
end;

function [h, wd] = textbox(x,y,str,c)
%  [h, wd] = textbox(x,y,str)    draws a box around the text 
% INPUT:  x, y - Coordinates
%         str  - Strings
% OUTPUT:    h - Object Handles
%           wd - x and y Width of boxes 

if ~isa(str,'cell'), str=cellstr(str); end
N = length(str);
wd = zeros(N,2);
h = zeros(N,2);
for i = 1:N,
    tx = text(x(i),y(i),str{i},'HorizontalAlignment','center','VerticalAlign','middle');
    sz = get(tx, 'Extent');
    wy = 2/3 * sz(4); wyB = y(i) - wy;  wyT = y(i) + wy;
    wx = max(2/3 * sz(3), wy); wxL = x(i) - wx; wxR = x(i) + wx;
    ptc = patch([wxL wxR wxR wxL], [wyT wyT wyB wyB], c(i,:)); 
    set(ptc, 'FaceColor', c(i,:));   %  'w' 
    wd(i,:) = [wx wy];
    delete(tx);
    tx = text(x(i),y(i),str{i},'HorizontalAlignment','center','VerticalAlign','middle');      
    h(i,:) = [tx ptc];
end;

function [h,yy,zz] = my_arrow(varargin)
% [h,yy,zz] = my_arrow(varargin)  Draw a line with an arrowhead.

% A lot of the original code is removed and most of the remaining can probably go too
% since it comes from a general use function only being called inone context. - Leon Peshkin 
% Copyright 1997, Erik A. Johnson <johnsone@uiuc.edu>, 8/14/97

ax         = [];       % set values to empty matrices
deflen        = 12;  %  16
defbaseangle  = 45;  %  90
deftipangle   = 16;
defwid = 0;  defpage = 0;  defends = 1;
ArrowTag = 'Arrow';  % The 'Tag' we'll put on our arrows
start      = varargin{1};    % fill empty arguments
stop       = varargin{2}; 
crossdir   = [NaN NaN NaN];   
len        = NaN; baseangle  = NaN;  tipangle = NaN;   wid = NaN;              
page       = 0; ends  = NaN;   
start = [start NaN];   stop = [stop NaN];
o         = 1;     % expand single-column arguments
ax        = gca;
% set up the UserData data (here so not corrupted by log10's and such)
ud = [start stop len baseangle tipangle wid page crossdir ends];
% Get axes limits, range, min; correct for aspect ratio and log scale
axm  = zeros(3,1);   axr = axm;   axrev = axm;  ap  = zeros(2,1);
xyzlog = axm; limmin    = ap;  limrange  = ap;  oldaxlims = zeros(1,7);
oneax = 1;      % all(ax==ax(1));  LPM
if (oneax),
	T = zeros(4,4); invT = zeros(4,4);
else
	T = zeros(16,1); invT = zeros(16,1); 
end
axnotdone = 1; % logical(ones(size(ax)));  LPM 
while (any(axnotdone))
	ii = 1;  % LPM min(find(axnotdone));
	curax = ax(ii);
	curpage = page(ii);
	% get axes limits and aspect ratio
	axl = [get(curax,'XLim'); get(curax,'YLim'); get(curax,'ZLim')];
	oldaxlims(find(oldaxlims(:,1)==0, 1),:) = [curax reshape(axl',1,6)];
	% get axes size in pixels (points)
	u = get(curax,'Units');
	axposoldunits = get(curax,'Position');
	really_curpage = curpage & strcmp(u,'normalized');
	if (really_curpage)
		curfig = get(curax,'Parent');  		pu = get(curfig,'PaperUnits');
		set(curfig,'PaperUnits','points');  pp = get(curfig,'PaperPosition');
		set(curfig,'PaperUnits',pu);         set(curax,'Units','pixels');
		curapscreen = get(curax,'Position'); set(curax,'Units','normalized');
		curap = pp.*get(curax,'Position');
    else
		set(curax,'Units','pixels');
		curapscreen = get(curax,'Position');
		curap = curapscreen;
    end
	set(curax,'Units',u);      set(curax,'Position',axposoldunits);
	% handle non-stretched axes position
	str_stretch = {'DataAspectRatioMode'; 'PlotBoxAspectRatioMode' ; 'CameraViewAngleMode' };
	str_camera  = {'CameraPositionMode'  ; 'CameraTargetMode' ; ...
	                'CameraViewAngleMode' ; 'CameraUpVectorMode'};
	notstretched = strcmp(get(curax,str_stretch),'manual');
	manualcamera = strcmp(get(curax,str_camera),'manual');
	if ~arrow_WarpToFill(notstretched,manualcamera,curax)
		% find the true pixel size of the actual axes
		texttmp = text(axl(1,[1 2 2 1 1 2 2 1]), ...
		               axl(2,[1 1 2 2 1 1 2 2]), axl(3,[1 1 1 1 2 2 2 2]),'');
		set(texttmp,'Units','points');
		textpos = get(texttmp,'Position');
		delete(texttmp);
		textpos = cat(1,textpos{:});
		textpos = max(textpos(:,1:2)) - min(textpos(:,1:2));
		% adjust the axes position
		if (really_curpage)  			% adjust to printed size
			textpos = textpos * min(curap(3:4)./textpos);
			curap = [curap(1:2)+(curap(3:4)-textpos)/2 textpos];
        else                         % adjust for pixel roundoff
			textpos = textpos * min(curapscreen(3:4)./textpos);
			curap = [curap(1:2)+(curap(3:4)-textpos)/2 textpos];
        end
    end
	% adjust limits for log scale on axes
	curxyzlog = [strcmp(get(curax,'XScale'),'log'); ...
	             strcmp(get(curax,'YScale'),'log'); strcmp(get(curax,'ZScale'),'log')];
	if (any(curxyzlog))
		ii = find([curxyzlog;curxyzlog]);
		if (any(axl(ii)<=0))
			error([upper(mfilename) ' does not support non-positive limits on log-scaled axes.']);
        else
			axl(ii) = log10(axl(ii));
        end
    end
	% correct for 'reverse' direction on axes;
	curreverse = [strcmp(get(curax,'XDir'),'reverse'); ...
	              strcmp(get(curax,'YDir'),'reverse'); strcmp(get(curax,'ZDir'),'reverse')];
	ii = find(curreverse);
	if ~isempty(ii)
		axl(ii,[1 2])=-axl(ii,[2 1]);
    end
	% compute the range of 2-D values
	curT = get(curax,'Xform');
	lim = curT*[0 1 0 1 0 1 0 1;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1;1 1 1 1 1 1 1 1];
	lim = lim(1:2,:)./([1;1]*lim(4,:));
	curlimmin = min(lim,[],2);
	curlimrange = max(lim,[],2) - curlimmin;
	curinvT = inv(curT);
	if ~oneax
        curT = curT.'; curinvT = curinvT.'; curT = curT(:); curinvT = curinvT(:);
    end
	% check which arrows to which cur corresponds
	ii = find((ax==curax)&(page==curpage));
	oo = ones(1,length(ii)); 	axr(:,ii) = diff(axl,1,2) * oo;
	axm(:,ii) = axl(:,1) * oo;  axrev(:,ii) = curreverse  * oo;
	ap(:,ii)  = curap(3:4)' * oo; xyzlog(:,ii) = curxyzlog   * oo;
	limmin(:,ii) = curlimmin  * oo;  limrange(:,ii) = curlimrange * oo;
	if (oneax),
		T    = curT;  invT = curinvT;
    else
		T(:,ii) = curT * oo; invT(:,ii) = curinvT * oo;
	end;
	axnotdone(ii) = zeros(1,length(ii));
end;
oldaxlims(oldaxlims(:,1)==0,:) = [];

% correct for log scales
curxyzlog = xyzlog.';  ii = find(curxyzlog(:));
if ~isempty(ii)
	start(ii) = real(log10(start(ii))); stop(ii) = real(log10(stop(ii)));
	if (all(imag(crossdir)==0)) % pulled (ii) subscript on crossdir, 12/5/96 eaj
		crossdir(ii) = real(log10(crossdir(ii)));
    end
end
ii = find(axrev.');    % correct for reverse directions
if ~isempty(ii)
	start(ii) = -start(ii);  stop(ii) = -stop(ii); crossdir(ii) = -crossdir(ii);
end
start  = start.';  stop  = stop.';   % transpose start/stop values
% take care of defaults, page was done above
ii = find(isnan(start(:)));  if ~isempty(ii),  start(ii) = axm(ii)+axr(ii)/2;  end;
ii = find(isnan(stop(:)));  if ~isempty(ii),  stop(ii) = axm(ii)+axr(ii)/2;  end;
ii = find(isnan(crossdir(:))); if ~isempty(ii),  crossdir(ii) = zeros(length(ii),1); end;
ii = find(isnan(len));  if ~isempty(ii),  len(ii) = ones(length(ii),1)*deflen; end;
baseangle(ii) = ones(length(ii),1)*defbaseangle;  tipangle(ii) = ones(length(ii),1)*deftipangle; 
wid(ii) = ones(length(ii),1) * defwid;   ends(ii) = ones(length(ii),1) * defends;
% transpose rest of values
len  = len.';  baseangle = baseangle.'; tipangle  = tipangle.'; wid = wid.';  
page = page.'; crossdir  = crossdir.';  ends = ends.'; ax   = ax.';

% for all points with start==stop, start=stop-(verysmallvalue)*(up-direction);
ii = find(all(start==stop));
if ~isempty(ii)
	% find an arrowdir vertical on screen and perpendicular to viewer
	%	transform to 2-D
		tmp1 = [(stop(:,ii)-axm(:,ii))./axr(:,ii);ones(1,length(ii))];
		if (oneax), twoD=T*tmp1;
        else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,ii).*tmp1;
		      tmp2=zeros(4,4*length(ii)); tmp2(:)=tmp1(:);
		      twoD=zeros(4,length(ii)); twoD(:)=sum(tmp2)'; 
        end
		twoD=twoD./(ones(4,1)*twoD(4,:));
	%	move the start point down just slightly
		tmp1 = twoD + [0;-1/1000;0;0]*(limrange(2,ii)./ap(2,ii));
	%	transform back to 3-D
		if (oneax), threeD=invT*tmp1;
        else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT(:,ii).*tmp1;
		      tmp2=zeros(4,4*length(ii)); tmp2(:)=tmp1(:);
		      threeD=zeros(4,length(ii)); threeD(:)=sum(tmp2)'; 
        end
		start(:,ii) = (threeD(1:3,:)./(ones(3,1)*threeD(4,:))).*axr(:,ii)+axm(:,ii);
end;
% compute along-arrow points
%	transform Start points
	tmp1 = [(start-axm)./axr; 1];
	if (oneax), X0=T*tmp1;
    else  tmp1 = [tmp1;tmp1;tmp1;tmp1]; tmp1=T.*tmp1;
	      tmp2 = zeros(4,4); tmp2(:)=tmp1(:);
	      X0=zeros(4,1); X0(:)=sum(tmp2)'; 
    end
	X0=X0./(ones(4,1)*X0(4,:));
%	transform Stop points
	tmp1=[(stop-axm)./axr; 1];
	if (oneax), Xf=T*tmp1;
    else  tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T.*tmp1;
	      tmp2=zeros(4,4); tmp2(:)=tmp1(:);
	      Xf=zeros(4,1); Xf(:)=sum(tmp2)'; 
    end
	Xf=Xf./(ones(4,1)*Xf(4,:));
%	compute pixel distance between points
	D = sqrt(sum(((Xf(1:2,:)-X0(1:2,:)).*(ap./limrange)).^2));
%	compute and modify along-arrow distances
	len1 = len;
	len2 = len - (len.*tan(tipangle/180*pi)-wid/2).*tan((90-baseangle)/180*pi);
	slen0 = 0; 	slen1 = len1 .* ((ends==2)|(ends==3));
	slen2 = len2 .* ((ends==2)|(ends==3));
	len0 = 0; len1  = len1 .* ((ends==1)|(ends==3));
	len2  = len2 .* ((ends==1)|(ends==3));
      ii = find((ends==1)&(D<len2));  	%	for no start arrowhead
	  if ~isempty(ii),
		  slen0(ii) = D(ii)-len2(ii);
	  end;
	  ii = find((ends==2)&(D<slen2));  	%	for no end arrowhead
	  if ~isempty(ii),
		  len0(ii) = D(ii)-slen2(ii);
	  end;
	len1  = len1  + len0;    len2 = len2  + len0;
	slen1 = slen1 + slen0; 	slen2 = slen2 + slen0;
 	% note:  the division by D below will probably not be accurate if both
 	%        of the following are true:
 	%           1. the ratio of the line length to the arrowhead
 	%              length is large
 	%           2. the view is highly perspective.
%	compute stoppoints
	tmp1 = X0.*(ones(4,1)*(len0./D))+Xf.*(ones(4,1)*(1-len0./D));
	if (oneax), tmp3 = invT*tmp1;
    else  tmp1 = [tmp1;tmp1;tmp1;tmp1]; tmp1 = invT.*tmp1;
	      tmp2 = zeros(4,4); tmp2(:) = tmp1(:);
	      tmp3 = zeros(4,1); tmp3(:) = sum(tmp2)'; 
    end
	stoppoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute tippoints
	tmp1=X0.*(ones(4,1)*(len1./D))+Xf.*(ones(4,1)*(1-len1./D));
	if (oneax), tmp3=invT*tmp1;
    else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,1); tmp3(:)=sum(tmp2)'; 
    end
	tippoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute basepoints
	tmp1=X0.*(ones(4,1)*(len2./D))+Xf.*(ones(4,1)*(1-len2./D));
	if (oneax), tmp3=invT*tmp1;
    else  tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,1); tmp3(:)=sum(tmp2)'; 
    end
	basepoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute startpoints
	tmp1=X0.*(ones(4,1)*(1-slen0./D))+Xf.*(ones(4,1)*(slen0./D));
	if (oneax), tmp3=invT*tmp1;
    else  tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4); tmp2(:) = tmp1(:);
	      tmp3=zeros(4,1); tmp3(:) = sum(tmp2)'; 
    end
	startpoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute stippoints
	tmp1=X0.*(ones(4,1)*(1-slen1./D))+Xf.*(ones(4,1)*(slen1./D));
	if (oneax), tmp3=invT*tmp1;
    else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1 = invT.*tmp1;
	      tmp2=zeros(4,4); tmp2(:)=tmp1(:); 
	      tmp3=zeros(4,1); tmp3(:)=sum(tmp2)'; 
    end
	stippoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute sbasepoints
	tmp1=X0.*(ones(4,1)*(1-slen2./D))+Xf.*(ones(4,1)*(slen2./D));
	if (oneax), tmp3=invT*tmp1;
    else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,1); tmp3(:)=sum(tmp2)'; 
    end
	sbasepoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;

% compute cross-arrow directions for arrows with NormalDir specified
if (any(imag(crossdir(:))~=0)),
	ii = find(any(imag(crossdir)~=0));
	crossdir(:,ii) = cross((stop(:,ii)-start(:,ii))./axr(:,ii), ...
	                       imag(crossdir(:,ii))).*axr(:,ii);
end;
basecross  = crossdir + basepoint;  % compute cross-arrow directions
tipcross   = crossdir + tippoint;  sbasecross = crossdir + sbasepoint;
stipcross  = crossdir + stippoint;
ii = find(all(crossdir==0)|any(isnan(crossdir)));
if ~isempty(ii),
	numii = length(ii);
	%	transform start points
		tmp1 = [basepoint(:,ii) tippoint(:,ii) sbasepoint(:,ii) stippoint(:,ii)];
		tmp1 = (tmp1-axm(:,[ii ii ii ii])) ./ axr(:,[ii ii ii ii]);
		tmp1 = [tmp1; ones(1,4*numii)];
		if (oneax), X0=T*tmp1;
        else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,[ii ii ii ii]).*tmp1;
		      tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
		      X0=zeros(4,4*numii); X0(:)=sum(tmp2)'; 
        end
		X0=X0./(ones(4,1)*X0(4,:));
	%	transform stop points
		tmp1 = [(2*stop(:,ii)-start(:,ii)-axm(:,ii))./axr(:,ii);ones(1,numii)];
		tmp1 = [tmp1 tmp1 tmp1 tmp1];
		if (oneax) Xf=T*tmp1;
        else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,[ii ii ii ii]).*tmp1;
		      tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
		      Xf=zeros(4,4*numii); Xf(:)=sum(tmp2)'; 
        end
		Xf=Xf./(ones(4,1)*Xf(4,:));
	%	compute perpendicular directions
		pixfact = ((limrange(1,ii)./limrange(2,ii)).*(ap(2,ii)./ap(1,ii))).^2;
		pixfact = [pixfact pixfact pixfact pixfact];
		pixfact = [pixfact;1./pixfact];
		[dummyval,jj] = max(abs(Xf(1:2,:)-X0(1:2,:)));
		jj1 = ((1:4)'*ones(1,length(jj))==ones(4,1)*jj);
		jj2 = ((1:4)'*ones(1,length(jj))==ones(4,1)*(3-jj));
		jj3 = jj1(1:2,:);
		Xp = X0;
		Xp(jj2) = X0(jj2) + ones(sum(jj2(:)),1);
		Xp(jj1) = X0(jj1) - (Xf(jj2)-X0(jj2))./(Xf(jj1)-X0(jj1)) .* pixfact(jj3);
	%	inverse transform the cross points
		if (oneax), Xp=invT*Xp;
		else, tmp1=[Xp;Xp;Xp;Xp]; tmp1=invT(:,[ii ii ii ii]).*tmp1;
		      tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
		      Xp=zeros(4,4*numii); Xp(:)=sum(tmp2)'; end;
		Xp=(Xp(1:3,:)./(ones(3,1)*Xp(4,:))).*axr(:,[ii ii ii ii])+axm(:,[ii ii ii ii]);
		basecross(:,ii)  = Xp(:,0*numii+(1:numii));
		tipcross(:,ii)   = Xp(:,1*numii+(1:numii));
		sbasecross(:,ii) = Xp(:,2*numii+(1:numii));
		stipcross(:,ii)  = Xp(:,3*numii+(1:numii));
end;

% compute all points
%	compute start points
	axm11 = [axm axm axm axm axm axm axm axm axm axm axm];
	axr11 = [axr axr axr axr axr axr axr axr axr axr axr];
	st = [stoppoint tippoint basepoint sbasepoint stippoint startpoint stippoint sbasepoint basepoint tippoint stoppoint];
	tmp1 = (st - axm11) ./ axr11;
	tmp1 = [tmp1; ones(1,size(tmp1,2))];
	if (oneax), X0=T*tmp1;
    else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[T T T T T T T T T T T].*tmp1;
	      tmp2=zeros(4,44); tmp2(:)=tmp1(:);
	      X0=zeros(4,11); X0(:)=sum(tmp2)'; 
    end
	X0=X0./(ones(4,1)*X0(4,:));
%	compute stop points
	tmp1 = ([start tipcross basecross sbasecross stipcross stop stipcross sbasecross basecross tipcross start] ...
	     - axm11) ./ axr11;
	tmp1 = [tmp1; ones(1,size(tmp1,2))];
	if (oneax), Xf=T*tmp1;
    else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[T T T T T T T T T T T].*tmp1;
	      tmp2=zeros(4,44); tmp2(:)=tmp1(:);
	      Xf=zeros(4,11); Xf(:)=sum(tmp2)'; 
    end
	Xf=Xf./(ones(4,1)*Xf(4,:));
%	compute lengths
	len0  = len.*((ends==1)|(ends==3)).*tan(tipangle/180*pi);
	slen0 = len.*((ends==2)|(ends==3)).*tan(tipangle/180*pi);
	le = [0 len0 wid/2 wid/2 slen0 0 -slen0 -wid/2 -wid/2 -len0 0];
	aprange = ap./limrange;
	aprange = [aprange aprange aprange aprange aprange aprange aprange aprange aprange aprange aprange];
	D = sqrt(sum(((Xf(1:2,:)-X0(1:2,:)).*aprange).^2));
	Dii=find(D==0); if ~isempty(Dii), D=D+(D==0); le(Dii)=zeros(1,length(Dii)); end; 
	tmp1 = X0.*(ones(4,1)*(1-le./D)) + Xf.*(ones(4,1)*(le./D));
%	inverse transform
	if (oneax), tmp3=invT*tmp1;
    else tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[invT invT invT invT invT invT invT invT invT invT invT].*tmp1;
	      tmp2=zeros(4,44); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,11); tmp3(:)=sum(tmp2)'; 
    end
	pts = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)) .* axr11 + axm11;
% correct for ones where the crossdir was specified
ii = find(~(all(crossdir==0)|any(isnan(crossdir))));
if ~isempty(ii),
	D1 = [pts(:,1+ii)-pts(:,9+ii) pts(:,2+ii)-pts(:,8+ii) ...
	      pts(:,3+ii)-pts(:,7+ii) pts(:,4+ii)-pts(:,6+ii) ...
	      pts(:,6+ii)-pts(:,4+ii) pts(:,7+ii)-pts(:,3+ii) ...
	      pts(:,8+ii)-pts(:,2+ii) pts(:,9+ii)-pts(:,1+ii)]/2;
	ii = ii'*ones(1,8) + ones(length(ii),1)*[1:4 6:9];   ii = ii(:)';
	pts(:,ii) = st(:,ii) + D1;
end;
% readjust for reverse directions
iicols = (1:1)';  iicols = iicols(:,ones(1,11));  iicols = iicols(:).';
tmp1 = axrev(:,iicols);
ii = find(tmp1(:)); if ~isempty(ii), pts(ii)=-pts(ii); end;
% readjust for log scale on axes
tmp1 = xyzlog(:,iicols);
ii = find(tmp1(:)); if ~isempty(ii), pts(ii)=10.^pts(ii); end;
% compute the x,y,z coordinates of the patches;
ii = (0:10)' + ones(11,1);
ii = ii(:)';
x = zeros(11,1);  y = x;    z = x;
x(:) = pts(1,ii)';   y(:) = pts(2,ii)';  z(:) = pts(3,ii)';
           % do the output
  % % create or modify the patches
H = 0; 
   % % make or modify the arrows
if arrow_is2DXY(ax(1)), zz=[]; else zz=z(:,1); end;
xyz = {'XData',x(:,1),'YData',y(:,1),'ZData',zz,'Tag',ArrowTag};
H(1) = patch(xyz{:});
  % % additional properties
set(H,'Clipping','off');
set(H,{'UserData'},num2cell(ud,2));
  % make sure the axis limits did not change

function [out,is2D] = arrow_is2DXY(ax)
% check if axes are 2-D X-Y plots,  may not work for modified camera angles, etc.
	out = zeros(size(ax)); % 2-D X-Y plots
	is2D = out;            % any 2-D plots
	views = get(ax(:),{'View'});
	views = cat(1,views{:});
	out(:) = abs(views(:,2))==90;
	is2D(:) = out(:) | all(rem(views',90)==0)';

function out = arrow_WarpToFill(notstretched,manualcamera,curax) %#ok<INUSL>
% check if we are in "WarpToFill" mode.
	out = strcmp(get(curax,'WarpToFill'),'on');
	% 'WarpToFill' is undocumented, so may need to replace this by
	% out = ~( any(notstretched) & any(manualcamera) );
