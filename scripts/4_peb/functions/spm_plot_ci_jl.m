function spm_plot_ci_jl(E,C,col, x,j,s)
% Plot mean and conditional confidence intervals
% FORMAT spm_plot_ci(E,C,x,j,s)
% E - expectation (structure or array)
% C - variance or covariance (structure or array)
% x - domain
% j - rows of E to plot
% s - string to specify plot type:e.g. '--r' or 'exp'
%
% If E is a row vector with two elements, confidence regions will be
% plotted; otherwise, bar charts with confidence intervals are provided
%__________________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_plot_ci.m 7305 2018-05-07 13:35:06Z karl $


%bar colour
%col=[0.8 0.8 0.8] %e.g. dark blue

% get axis
%--------------------------------------------------------------------------
ax = gca;

% confidence region (CR) plotting
%--------------------------------------------------------------------------
if size(E,1) == 1 && size(E,2) == 2
    E  = E';
    CR = true;
else
    CR = false;
end
    

% unpack expectations into a matrix
%--------------------------------------------------------------------------
if isstruct(E),       E = spm_vec(E);    end
if iscell(E),         E = spm_vec(E);    end

if ~exist('x','var'), x = 1:size(E,2);   end
if ~exist('j','var'), j = 1:size(E,1);   end
if ~exist('s','var'), s = '';            end

if isempty(x),        x = 1:size(E,2);   end
if isempty(j),        j = 1:size(E,1);   end


% order and length of sequence
%--------------------------------------------------------------------------
O     = E;
E     = E(j,:);
[n,N] = size(E);

% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);               % confidence interval
gr    = 0.9;                                 % grey level
if iscell(C)
    
    % try cell array of covariances (from spm_DEM amd spm_LAP)
    %----------------------------------------------------------------------
    try
        for i = 1:N
            c(:,i) = ci*sqrt(diag(C{i}(j,j)));
        end
    catch
        
        % try cell array of variances
        %------------------------------------------------------------------
        c = ci*sqrt(spm_unvec(spm_vec(C),O));
        c = c(j,:);
    end
    
elseif isstruct(C)
    
    % try structure of variances
    %----------------------------------------------------------------------
    c = ci*sqrt(spm_unvec(spm_vec(C),O));
    c = c(j,:);
    
elseif isnumeric(C)
    
    % try matrix of variances
    %----------------------------------------------------------------------
    if all(size(C) == size(O))
        c = ci*sqrt(C(j,:));
    elseif all(size(C') == size(O))
        c = ci*sqrt(C(:,j));
    else
        
        % try covariance matrix
        %------------------------------------------------------------------
        C = diag(C);
        c = ci*sqrt(C(j,:));
    end
    
end

% set plot parameters
%--------------------------------------------------------------------------
switch lower(get(ax,'NextPlot'))
    case 'add'
        
        
        width = .9;
    otherwise
        
        width = .8;
        
end

% plot elliptical confidence region
%--------------------------------------------------------------------------
if CR
    [x,y] = ellipsoid(E(1),E(2),1,c(1),c(2),0,32);
    fill(x(16,:)',y(16,:)',[1 1 1]*gr,'EdgeColor',[1 1 1]*.5,'Parent',ax);
    hold(ax,'on');
    plot(ax,E(1),E(2),'.','MarkerSize',16);
    hold(ax,'off'); drawnow
    return
end


% plot bar chart
%--------------------------------------------------------------------------
if N >= 8
    
    % time-series plot
    %======================================================================
    if strcmpi(s,'exp')
        fill([x fliplr(x)],exp([full(E + c) fliplr(full(E - c))]),...
            [.95 .95 1],'EdgeColor',[.8 .8 1],'Parent',ax);
        hold(ax,'on');
        plot(x,exp(E));
        
    else
        fill([x fliplr(x)],[full(E + c) fliplr(full(E - c))],...
            [.95 .95 1],'EdgeColor',[.8 .8 1],'Parent',ax);
        hold(ax,'on');
        plot(ax,x,E,s);
    end    
    
else
    
    
    % bar
    %======================================================================
    if N == 1
        
        if strcmpi(s,'exp')
            
            % conditional means
            %--------------------------------------------------------------
            bar(ax,exp(E),width,'Edgecolor',[1 1 1]/2,'Facecolor',[0 0.4471 0.7412]*.8);
            hold(ax,'on');
            
            % conditional variances
            %--------------------------------------------------------------
            for k = 1:n
                line([k k],exp([-1 1]*c(k) + E(k)),'LineWidth',4,'Color',col,'Parent',ax);
            end
            
            
            
        else
                        
            if n > 1
                
                % conditional means
                %----------------------------------------------------------
                bar(ax,E,width,'Edgecolor',[1 1 1]/2,'Facecolor',col);
                hold(ax,'on');
                
            else
                % conditional means
                %----------------------------------------------------------
                bar(ax,E,'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8, 'LineWidth',2);
                hold(ax,'on');
                
            end
            
            % conditional variances
            %--------------------------------------------------------------
            for k = 1:n
                line([k k],[-1 1]*c(k) + E(k),'LineWidth',2,'Color',[0.8 0.8 0.8],'Parent',ax);
            end
            
        end
        
        box(ax,'off');
        set(ax,'XLim',[0 n + 1]);
        
    else
        
        if strcmpi(s,'exp')
            
            % conditional means (exponential)
            %------------------------------------------------------------------
            h = bar(ax,exp(E)'); hold(ax,'on');
            
            % conditional variances
            %------------------------------------------------------------------
            for m = 1:n
                x = mean(get(get(h(m),'Children'),'Xdata'));
                for k = 1:N
                    line([x(k) x(k)],exp([-1 1]*c(m,k) + E(m,k)),'LineWidth',1,'Color',col,'Parent',ax);
                end
            end
            
        else
            
            % conditional means
            %------------------------------------------------------------------
            h = bar(ax,E); hold(ax,'on');
            
            % conditional variances
            %------------------------------------------------------------------
            for m = 1:N
                x = mean(get(get(h(m),'Children'),'Xdata'));
                for k = 1:n
                    line([x(k) x(k)],[-1 1]*c(k,m) + E(k,m),'LineWidth',4,'Color',col,'Parent',ax);
                end
            end
            
        end
    end
    
end
hold(ax,'off');
drawnow
