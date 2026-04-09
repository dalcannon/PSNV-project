function demo( X, wave )
%---------------------------------------------------------------------
% 
% PURPOSE: Demo performing all spectral transformations. 
% See [A] and [B] below. 
% 
%---------------------------------------------------------------------
% 
% USAGE
% regional_demo( X, wave )
% 
%---------------------------------------------------------------------
% 
% INPUT
% [1] X: (m,n) matrix of spectra.  Each spectrum is aligned row-wise.
%
% [2] wave: n-vector. (OPTIONAL.)  Wavelengths.
%
%---------------------------------------------------------------------
%
% [A] MULTIPLICATIVE SCATTER CORRECTION AND VARIANTS
% 
% [A1] msc.m: Multiplicative Scatter Correction (MSC)
% 
% [A2] lmsc.m: Local MSC
% 
% [A3] pmsc.m: Moving Piecewise MSC (no loops)
%      Only implements Helland version of MSC.
% 
%---------------------------------------------------------------------
%
% [B] STANDARD NORMAL VARIATE CORRECTION AND VARIANTS
% 
% [B1] snv.m: Standard Normal Variate (SNV) correction
% 
% [B2] lsnv.m: Local SNV
% 
% [B4] psnv.m: Moving Piecewise SNV (no loops)
%
%---------------------------------------------------------------------

%-----------------------------------------------------------
% Sanity checks on spectra 
%-----------------------------------------------------------
% Input validation and preprocessing
if ~ismatrix(X)
    error('X must be a matrix');
end
% Convert to matrix if vector
if isvector(X)
    X = X(:)'; 
end
[~, n] = size(X); 


%-----------------------------------------------------------
% Sanity checks on wavelengths
%-----------------------------------------------------------
if nargin < 2
    wave = 1:n;
end   
wave = sort(wave(:))';
if length(wave) ~=n
    error('wave must have %d entries',n);
end    


%-----------------------------------------------------------
% Default parameters
% w: Window size used in piecewise methods
% r: Reference spectrum
% k: Number of disjoint segments/blocks used in local methods
%-----------------------------------------------------------
w = 31;
r = mean(X,1);  
k = 25;


%-----------------------------------------------------------
% Set y-axis limits for each subplot 
%-----------------------------------------------------------
yrange = {[],[],[],[]};


%-----------------------------------------------------------
% Original spectra
%-----------------------------------------------------------
% 
% Plot spectra
%
axes('Position',[0.06 0.55 0.41 0.41]);
plot(wave,X,'k');
gridlines(wave);
if ~isempty(yrange{1}), set(gca,'YLim',yrange{1}); end
title('ORIGINAL SPECTRA','FontSize',18);

%-----------------------------------------------------------
% Helland-based MSC and variants
%-----------------------------------------------------------
%
% Transform spectra 
%
classic = false;
Z = struct();
Z.msc = msc( X, r, classic );
Z.lmsc = lmsc( X, k, r, classic );
Z.mpmsc = pmsc( X, w, r );
%
% Plot spectra
%
axes('Position',[0.06 0.06 0.41 0.41]);
plot_variants(Z,wave);
xlabel('WAVELENGTH','FontSize',16);
if ~isempty(yrange{3}), set(gca,'YLim',yrange{3}); end
title('HELLAND MSC AND VARIANTS','FontSize',18);


%-----------------------------------------------------------
% SNV and variants
%-----------------------------------------------------------  
%
% Transform spectra 
%
Z = struct();
Z.snv = snv( X );
Z.lsnv = lsnv( X, k );
Z.psnv = psnv( X, w );
%
% Plot spectra
%
axes('Position',[0.53 0.06 0.41 0.41]);
plot_variants(Z,wave);
xlabel('WAVELENGTH','FontSize',16);
if ~isempty(yrange{4}), set(gca,'YLim',yrange{4}); end
title('CLASSICAL SNV AND VARIANTS','FontSize',18);


%---------------------------------------------------------------------
end


function plot_variants(Xtransform,wave)
%---------------------------------------------------------------------
F = fieldnames(Xtransform);
N = length(F);
Zcolor = jet(N);
for i = 1:N
    fprintf('METHOD = %s\n',F{i});
    Z = Xtransform.(F{i});
    h = plot(wave,Z);
    set(h,'color',Zcolor(i,:),'linewidth',2);
    if i==1
        hlegend = repmat(h(1),N,1);
        hold on;
    else     
        hlegend(i) = h(1);
    end   
end    
hold off;
method = replace(upper(F), '_', '-');
legend(hlegend,method,'Location','southeast');
gridlines(wave);
%---------------------------------------------------------------------
end

function gridlines(wave)
%---------------------------------------------------------------------
set(gca,'XLim',[wave(1) wave(end)]);
set(gca,'ticklength',[0 0]);
set(gca,'color',0.55*[1 1 1]);
grid on;
set(gca,'GridColor',0.95*[1 1 1]);
set(gca,'GridAlpha',0.90);
set(gca,'GridLineStyle','-');
%---------------------------------------------------------------------
end
