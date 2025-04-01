function demo(dataset)
%---------------------------------------------------------------------
% [S]pectral [T]ransformation [DEMO]: Perform and compare various 
% MSC and SNV variants.  Each method has its own acronym (and 
% corresponding m-file) in single quotes.
% 
% MSC-BASED TRANSFORMATIONS
% 'msc':   MULTIPLICATIVE SCATTERT CORRECTION
% 'pmsc':  PIECEWISE MSC
% 'mpmsc': MOVING PIECEWISE MSC
% 'emsc':  EXTENDED MSC
% 'pemsc': PIECEWISE EXTENDED MSC
% 'lmsc':  LOCAL MSC
% 'loopy': LOOPY MSC
% 
% SNV-BASED TRANSFORMATIONS
% 'snv':   STANDARD NORMAL VARIATE
% 'psnv':  PIECEWISE SNV
% 'mpsnv': MOVING PIECEWISE SNV
% 'lsnv':  LOCAL SNV
% 
%---------------------------------------------------------------------
% 
% INPUT: 
% [1] dataset: string variable indiating which data set to load.  
%     The data sets currently include 'corn' and 'inor'.
%     From this data set, we will extract/create the variables: 
%     > 'wave': p-vector of continuous features such as wavelengths.
%     > 'X': (n,p) matrix of spectra (n samples by p features).
% 
%---------------------------------------------------------------------

%-------------------------------------------------
% LOAD DATA SET: GET WAVELENGTHS AND SPECTRA
%-------------------------------------------------
if nargin<1, dataset = 'corn'; end
switch dataset
    case 'corn'
        load('corn.mat','corn');
        wave = corn.wave;
        X = corn.Xm5;
    case 'inor'
        load('inor.mat','inor');
        wave = inor.wave;
        X = inor.X;      
    %case 'yourdataset'    
        %
    otherwise
        error('Invaldid data set');
end    
% Check length of wave
p = size(X,2);
wave = wave(:);
p2 = length(wave); 
if p2 ~= p
    error('length(wave) does not equal size(X,2)'); 
end 

%-------------------------------------------------
% MSC AND SNV DEMOS
%-------------------------------------------------
msc_demo(wave,X);
snv_demo(wave,X);
   
%---------------------------------------------------------------------
end

function msc_demo(wave,X)
%---------------------------------------------------------------------

p = size(X,2);
r = mean(X,1);
v = 25; % half-window for piecewise methods

%-------------------------------------------------
% We will compute MSC in two different 
% but equivalent ways
%-------------------------------------------------
% - Classical MSC
t0 = tic; Z1 = msc( X, r, true  ); t1 = toc(t0);
% - Helland-based MSC (see msc.m)
t0 = tic; Z2 = msc( X, r, false ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('MSC (%1.5f)',t1);
str2 = sprintf('Helland MSC (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});
% - Display log10(Frobenius norm) between Z1 and Z2
dZ = Z1-Z2;
err = log10(norm(dZ(:))); 
title(sprintf('log_{10}(||Helland-MSC - MSC||_F) = %2.4f',err));


%-------------------------------------------------
% We will compute PMSC in two different but 
% equivalent ways but the second way is faster
%-------------------------------------------------
% - Standard PMSC
t0 = tic; Z1 = pmsc(  X, v, r, true ); t1 = toc(t0);
% - Moving PMSC
t0 = tic; Z2 = mpmsc( X, v, r       ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('PMSC (%1.5f)',t1);
str2 = sprintf('MOVING PMSC (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});
% - Display log10(Frobenius norm) between Z1 and Z2
dZ = Z1-Z2;
err = log10(norm(dZ(:))); 
title(sprintf('log_{10}(||Helland-MSC - MSC||_F) = %2.4f',err));

%-------------------------------------------------
% We will compute EMSC in two different ways.  
% The first is the classical approach, while 
% the second is projection-based.
%-------------------------------------------------
k = 3; % typical type of emsc expansion (see emsc.m)
% - Classical EMSC
t0 = tic; Z1 = emsc( X, r, k, true );  t1 = toc(t0);
% - Projection-based EMSC (mean and linear trend removed)
t0 = tic; Z2 = emsc( X, r, k, false ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('EMSC (%1.5f)',t1);
str2 = sprintf('Projection-based PMSC (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});

%-------------------------------------------------
% We will compare EMSC versus Piecewise EMSC
%-------------------------------------------------
k = 2; % truncated type of emsc expansion (see emsc.m)
% - Classical EMSC
t0 = tic; Z1 = emsc(  X,    r, k, true ); t1 = toc(t0);
% - Piecewise EMSC (or PEMSC)
t0 = tic; Z2 = pemsc( X, v, r, k, true ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('EMSC (%1.5f)',t1);
str2 = sprintf('PIECEWISE PMSC (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});

%-------------------------------------------------
% MSC vs Loopy MSC with 8 iterations
%-------------------------------------------------
% - Classical MSC
t0 = tic; Z1 = msc( X, r, true  ); t1 = toc(t0);
% - Loopy MSC
N = 8;
t0 = tic; [Z2,RRSSQ] = loopy( X, v, N, 'msc' ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('MSC (%1.5f)',t1);
str2 = sprintf('LOOPY MSC (N=%d, %1.5f)',N,t2);
specplot2(wave,Z1,Z2{end},{str1,str2});
% 
figure; 
semilogy(1:N,RRSSQ,'-o');
xlabel('Loopy iterations'); ylabel('RRSSQ');
grid on; 

%-------------------------------------------------
% Compare MSC vs Local MSC
%-------------------------------------------------
% - Classical MSC
t0 = tic; Z1 = msc( X, r, true  ); t1 = toc(t0);
% - Local MSC
N = 10; % half-window (v) below will be ignored
t0 = tic; Z2 = lmsc( X, N, v, r, 'msc' ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('MSC (%1.5f)',t1);
str2 = sprintf('Local MSC (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});
itick = [1 round(p/N:p/N:p)];
set(gca,'xtick',wave(itick));


%---------------------------------------------------------------------
end

function snv_demo(wave,X)
%---------------------------------------------------------------------

p = size(X,2);
v = 25; % half-window for piecewise methods

%-------------------------------------------------
% Compare SNV vs PSNV
%-------------------------------------------------
% - SNV
t0 = tic; Z1 = snv(  X    ); t1 = toc(t0);
% - Piecewise SNV
t0 = tic; Z2 = psnv( X, v ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('SNV (%1.5f)',t1);
str2 = sprintf('Piecewise SNV (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});
set(gca,'color',[1 1 0.9]);


%-------------------------------------------------
% Compare PSNV vs moving PSNV
%-------------------------------------------------
% - Standard PSNV
t0 = tic; Z1 = psnv(  X, v ); t1 = toc(t0);
% - Moving PSNV
t0 = tic; Z2 = mpsnv( X, v ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('PSNV (%1.5f)',t1);
str2 = sprintf('MOVING PSNV (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2});
% - Display log10(Frobenius norm) between Z1 and Z2
dZ = Z1-Z2;
err = log10(norm(dZ(:))); 
title(sprintf('log_{10}(||Helland-MSC - MSC||_F) = %2.4f',err));
set(gca,'color',[1 1 0.9]);

%-------------------------------------------------
% Compare SNV vs Local SNV
%-------------------------------------------------
% - SNV
t0 = tic; Z1 = snv(  X    ); t1 = toc(t0);
% - Local SNV
N = 10; % Number of windows
t0 = tic; Z2 = lsnv( X, N ); t2 = toc(t0);
% - Plot corrected spectra and elapsed times
str1 = sprintf('SNV (%1.5f)',t1);
str2 = sprintf('Local SNV (%1.5f)',t2);
specplot2(wave,Z1,Z2,{str1,str2},false);
itick = [1 round(p/N:p/N:p)];
set(gca,'xtick',wave(itick));
set(gca,'color',[1 1 0.9]);

%---------------------------------------------------------------------
end

function specplot2(wave,Z1,Z2,str,plotline)
%---------------------------------------------------------------------
if (nargin<5), plotline = true; end
figure;
if plotline
    % Plot lines
    h = plot(wave,Z1,'b-',wave,Z2,'r-');
else
    % Plot dots instead
    h = plot(wave,Z1,'b.',wave,Z2,'r.');
end    
% Plot the x-axis
ylim = get(gca,'YLim');
if ylim(1)==0
    dy = diff(ylim)/20;
    set(gca,'ylim',[-dy ylim(2)]);
end    
hold on;
xlim = get(gca,'xlim');
plot(xlim,[0 0],'k-','linewidth',1);
hold off;
H = [h(1); h(end)];
grid;
legend(H,str,'Location','northwest');
%---------------------------------------------------------------------
end
