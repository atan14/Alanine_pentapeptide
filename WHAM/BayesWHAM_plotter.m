function BayesWHAM_plotter(f__hist_binCenters, f__hist_binWidths, f__pdf_MAP, f__betaF_MAP, f__f_MAP, f__pdf_MH, f__betaF_MH, f__f_MH, f__logL_MH, f__step_MH)

% Copyright:	Andrew L. Ferguson, UIUC 
% Last updated:	2 Jan 2016

% SYNOPSIS
%
% code to perform plotting of: 	(i)   {1,2,3}-dimensional maximum a posteriori (MAP) free energy surfaces (FES) computed by Bayesian inference of biased umbrella sampling trajectories in umbrella variables psi 
%								(ii)  uncertainty estimates in MAP FES estimated by Metropolis-Hastings (MH) sampling of Bayes posterior distribution 
% 								(iii) likelihood trajectory along MH sampling path to assess convergence of MH sampling 

% INPUTS
%
% f__hist_binCenters  		- [str] [dim x M_k float] path to text file containing dim rows specifying M_k k=1..dim centers of the rectilinear histogram bins in each dimension used to construct dim-dimensional histograms constituting f__pdf_MAP/MH and f__betaF_MAP/MH 
% f__hist_binWidths   		- [str] [dim x M_k float] path to text file containing dim rows specifying M_k k=1..dim widths of the rectilinear histogram bins in each dimension used to construct dim-dimensional histograms constituting f__pdf_MAP/MH and f__betaF_MAP/MH 
% f__pdf_MAP          		- [str] [1 x M float] path to text file containing MAP estimate of unbiased probability density function pdf_l over l=1..M bins of dim-dimensional rectilinear histogram in row major order (last index changes fastest) 
% f__betaF_MAP        		- [str] [1 x M float] path to text file containing MAP estimate of unbiased free energy surface beta*F_l = -ln(p(psi)/binVolume) + const. over l=1..M bins of dim-dimensional rectilinear histogram in row major order (last index changes fastest) 
% f__f_MAP            		- [str] [1 x S float] path to text file containing MAP estimates of f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i, for i=1..S biased simulations 
% f__pdf_MH           		- [str] [nSamples_MH x M float] path to text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of unbiased probability density function pdf_l over l=1..M bins of dim-dimensional rectilinear histogram in row major order (last index changes fastest) 
% f__betaF_MH         		- [str] [nSamples_MH x M float] path to text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of unbiased free energy surface beta*F_l = -ln(p(psi)/binVolume) + const. over l=1..M bins of dim-dimensional rectilinear histogram in row major order (last index changes fastest) 
% f__f_MH             		- [str] [nSamples_MH x S float] path to text file containing nSamples_MH Metropolis-Hastings samples from the Bayes posterior of f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i, for i=1..S biased simulations 
% f__logL_MH          		- [str] [nSamples_MH x 1 float] path to text file containing log of Bayes posterior up to an additive constant over the nSamples_MH Metropolis-Hastings samples from the posterior 
% f__step_MH          		- [str] [nSamples_MH x 1 int] path to text file containing Metropolis-Hastings step associated with each of the nSamples_MH Metropolis-Hastings samples from the posterior 

% OUTPUTS
%
% logL.fig/jpg				- trajectory of the log Bayes posterior up to an additive constant over the nSamples_MH Metropolis-Hastings samples from the Bayesian posterior 
% f.fig/jpg					- MAP estimates of f_i = Z/Z_i = ratio of unbiased partition function to that of biased simulation i for i=1..S biased simulations, overlaid with traces of MH samples of f_i 
%
% 1-dimension:
%
% pdf.fig/jpg				- MAP estimate of probability density function as a function of psi 
% pdf_traces.fig/jpg		- MAP estimate of probability density function as a function of psi overlaid with traces of MH samples 
% pdf_eb_limits.fig/jpg		- MAP estimate of probability density function as a function of psi overlaid with the MAP estimates plus and minus the standard deviation of the ensemble of MH samples 
% pdf_eb_stdev.fig/jpg		- standard deviation of the probability density function of the ensemble of MH samples as a function of psi 
% pdf_eb_bars.fig/jpg		- MAP estimate of probability density function as a function of psi with error bars denoting the standard deviation of the ensemble of MH samples 
% betaF.fig/jpg				- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi 
% betaF_traces.fig/jpg		- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi overlaid with traces of MH samples 
% betaF_eb_limits.fig/jpg	- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi overlaid with the MAP estimates plus and minus the standard deviation of the ensemble of MH samples 
% betaF_eb_stdev.fig/jpg	- standard deviation of the beta*F = -ln(p(psi)/binVolume) + const. of the ensemble of MH samples as a function of psi 
% betaF_eb_bars.fig/jpg		- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi with error bars denoting the standard deviation of the ensemble of MH samples 
% betaF_eb_dist.fig/jpg		- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi overlaid with the probability density function of the MH ensemble at each MAP data point 
%
% 2-dimensions:
%
% pdf.fig/jpg				- MAP estimate of probability density function as a function of psi 
% pdf_traces.fig/jpg		- MAP estimate of probability density function as a function of psi overlaid with traces of MH samples 
% pdf_eb_limits.fig/jpg		- MAP estimate of probability density function as a function of psi overlaid with the MAP estimates plus and minus the standard deviation of the ensemble of MH samples 
% pdf_eb_stdev.fig/jpg		- standard deviation of the probability density function of the ensemble of MH samples as a function of psi 
% betaF.fig/jpg				- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi 
% betaF_traces.fig/jpg		- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi overlaid with traces of MH samples 
% betaF_eb_limits.fig/jpg	- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi overlaid with the MAP estimates plus and minus the standard deviation of the ensemble of MH samples 
% betaF_eb_stdev.fig/jpg	- standard deviation of the beta*F = -ln(p(psi)/binVolume) + const. of the ensemble of MH samples as a function of psi 
%
% 3-dimensions:
%
% pdf.fig/jpg				- MAP estimate of probability density function as a function of psi 
% pdf_eb_stdev.fig/jpg		- standard deviation of the probability density function of the ensemble of MH samples as a function of psi 
% betaF.fig/jpg				- MAP estimate of beta*F = -ln(p(psi)/binVolume) + const. as a function of psi 
% betaF_eb_stdev.fig/jpg	- standard deviation of the beta*F = -ln(p(psi)/binVolume) + const. of the ensemble of MH samples as a function of psi 


% parameters
fsize=16;


% default args
if nargin == 0
    f__hist_binCenters = 'hist_binCenters.txt';
    f__hist_binWidths = 'hist_binWidths.txt';
    f__pdf_MAP = 'pdf_MAP.txt';
    f__betaF_MAP = 'betaF_MAP.txt';
    f__f_MAP = 'f_MAP.txt';
    f__pdf_MH = 'pdf_MH.txt';
    f__betaF_MH = 'betaF_MH.txt';
    f__f_MH = 'f_MH.txt';
    f__logL_MH = 'logL_MH.txt';
    f__step_MH = 'step_MH.txt';
end


% loading data
fprintf('\n');
fprintf('Loading data...\n');

binC = cell(0);
fin = fopen(f__hist_binCenters,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        binC{ii} = transpose(C{1});
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

binW = cell(0);
fin = fopen(f__hist_binWidths,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        binW{ii} = transpose(C{1});
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

pdf_MAP = [];
fin = fopen(f__pdf_MAP,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        pdf_MAP(ii,:) = C{1};
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

betaF_MAP = [];
fin = fopen(f__betaF_MAP,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        betaF_MAP(ii,:) = C{1};
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

f_MAP = [];
fin = fopen(f__f_MAP,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        f_MAP(ii,:) = C{1};
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

pdf_MH = [];
fin = fopen(f__pdf_MH,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        pdf_MH(ii,:) = C{1};
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

betaF_MH = [];
fin = fopen(f__betaF_MH,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        betaF_MH(ii,:) = C{1};
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

f_MH = [];
fin = fopen(f__f_MH,'rt');
    ii=1;
	tline = fgetl(fin);
    while ischar(tline)
        C = textscan(tline,'%f');
        f_MH(ii,:) = C{1};
        ii=ii+1;
        tline = fgetl(fin);
    end
fclose(fin);

fin = fopen(f__logL_MH,'rt');
    C = textscan(fin,'%f');
    logL_MH = C{1};
fclose(fin);

fin = fopen(f__step_MH,'rt');
    C = textscan(fin,'%f');
    step_MH = C{1};
fclose(fin);

fprintf('DONE!\n\n')


% post-processing and error checking
fprintf('Error checking data import...\n')

dim = length(binC);
if length(binW) ~= dim
    error('Dimensionality of %s and %s are incompatible',f__hist_binCenters,f__hist_binWidths);
end
    
M_k = nan(dim,1);
for ii=1:dim
    if size(binC{ii},1) ~= 1
        error('Dimension %d of %s is not a row vector',ii,f__hist_binCenters);
    end
    if size(binW{ii},1) ~= 1
        error('Dimension %d of %s is not a row vector',ii,f__hist_binWidths);
    end
    M_k(ii) = length(binC{ii});
    if length(binW{ii}) ~= M_k(ii)
        error('Number of bins in dimension %d of %s and %s are incompatible',ii,f__hist_binCenters,f__hist_binWidths);
    end
end

M = prod(M_k);
if size(pdf_MAP,1) ~= 1
    error('File %s does not contain a row vector',f__pdf_MAP);
end
if size(pdf_MAP,2) ~= M
    error('Number of bins in %s is incompatible with %s and %s',f__pdf_MAP,f__hist_binCenters,f__hist_binWidths);
end
if size(betaF_MAP,1) ~= 1
    error('File %s does not contain a row vector',f__betaF_MAP);
end
if size(betaF_MAP,2) ~= M
    error('Number of bins in %s is incompatible with %s and %s',f__betaF_MAP,f__hist_binCenters,f__hist_binWidths);
end

S = size(f_MAP,2);
if size(f_MAP,1) ~= 1
    error('File %s does not contain a row vector',f__f_MAP);
end

nSamples_MH = size(pdf_MH,1);
if size(pdf_MH,2) ~= M
    error('Number of bins in %s is incompatible with %s and %s',f__pdf_MH,f__hist_binCenters,f__hist_binWidths);
end
if size(betaF_MH,1) ~= nSamples_MH
    error('Number of rows (i.e., M-H samples) in %s is incommensurate with that in %s',f__betaF_MH,f__pdf_MH);
end
if size(betaF_MH,2) ~= M
    error('Number of bins in %s is incompatible with %s and %s',f__betaF_MH,f__hist_binCenters,f__hist_binWidths);
end
if size(f_MH,1) ~= nSamples_MH
    error('Number of rows (i.e., M-H samples) in %s is incommensurate with that in %s',f__f_MH,f__pdf_MH);
end
if size(f_MH,2) ~= S
    error('Number of columns (i.e., simulations) in %s is incompatible with %s',f__f_MH,f__f_MAP);
end
if size(logL_MH,1) ~= nSamples_MH
    error('Number of rows (i.e., M-H samples) in %s is incommensurate with that in %s',f__logL_MH,f__pdf_MH);
end
if size(logL_MH,2) ~= 1
    error('File %s does not contain a column vector',f__logL_MH);
end
if size(step_MH,1) ~= nSamples_MH
    error('Number of rows (i.e., M-H samples) in %s is incommensurate with that in %s',f__step_MH,f__pdf_MH);
end
if size(step_MH,2) ~= 1
    error('File %s does not contain a column vector',f__step_MH);
end

fprintf('DONE!\n\n')


% plotting
fprintf('Plotting...\n');

% log likelihood
figure;
plot(step_MH,logL_MH,'k');
xlabel('M-H step','fontsize',fsize);
ylabel('log(L) / -','fontsize',fsize);
set(gca,'fontsize',fsize);
set(gcf,'color','w');
saveas(gcf,'logL','fig');
saveas(gcf,'logL','jpg');

% f
figure;
semilogy( transpose( repmat(1:length(f_MAP),size(f_MH,1),1) ), transpose(f_MH), 'color',[0.5 0.5 0.5] );
hold on
semilogy(1:S,f_MAP,'r');
hold off
xlabel('biased simulation i','fontsize',fsize);
ylabel('f_i / -','fontsize',fsize);
set(gca,'fontsize',fsize);
set(gcf,'color','w');
saveas(gcf,'f','fig');
saveas(gcf,'f','jpg');

if dim==1
    
    % pdf
    
    % - naked
    figure;
    plot(binC{1},pdf_MAP,'-r');
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('pdf(\psi_1) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf','fig');
    saveas(gcf,'pdf','jpg');
    
    % - traces
    figure;
    plot( transpose( repmat(binC{1},size(pdf_MH,1),1) ), transpose(pdf_MH), 'color',[0.5 0.5 0.5] );
    hold on
    plot(binC{1},pdf_MAP,'-r');
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('pdf(\psi_1) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_traces','fig');
    saveas(gcf,'pdf_traces','jpg');
    
    % - stdev errorbars
    pdf_MH_std = std(pdf_MH,0,1);
    
    figure;
    plot(binC{1},pdf_MAP,'-r');
    hold on
    plot(binC{1},pdf_MAP+pdf_MH_std,'color',[0.5 0.5 0.5]);
    plot(binC{1},pdf_MAP-pdf_MH_std,'color',[0.5 0.5 0.5]);
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('pdf(\psi_1) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_eb_limits','fig');
    saveas(gcf,'pdf_eb_limits','jpg');
    
    figure;
    plot(binC{1},pdf_MH_std,'-b');
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('stdev(pdf(\psi_1)) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_eb_stdev','fig');
    saveas(gcf,'pdf_eb_stdev','jpg');
    
    figure;
    errorbar(binC{1},pdf_MAP,pdf_MH_std,'-r');
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('pdf(\psi_1) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_eb_bars','fig');
    saveas(gcf,'pdf_eb_bars','jpg');
    
    
    % betaF
    
    % - naked
    figure;
    plot(binC{1},betaF_MAP,'-r');
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\betaF(\psi_1) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF','fig');
    saveas(gcf,'betaF','jpg');
    
    % - traces
    figure;
    plot( transpose( repmat(binC{1},size(betaF_MH,1),1) ), transpose(betaF_MH), 'color',[0.5 0.5 0.5] );
    hold on
    plot(binC{1},betaF_MAP,'-r');
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\betaF(\psi_1) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_traces','fig');
    saveas(gcf,'betaF_traces','jpg');
    
    % - stdev errorbars
    betaF_MH_std = std(betaF_MH,0,1);
    
    figure;
    plot(binC{1},betaF_MAP,'-r');
    hold on
    plot(binC{1},betaF_MAP+betaF_MH_std,'color',[0.5 0.5 0.5]);
    plot(binC{1},betaF_MAP-betaF_MH_std,'color',[0.5 0.5 0.5]);
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\betaF(\psi_1) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_limits','fig');
    saveas(gcf,'betaF_eb_limits','jpg');
    
    figure;
    plot(binC{1},betaF_MH_std,'-b');
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('stdev(\betaF(\psi_1)) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_stdev','fig');
    saveas(gcf,'betaF_eb_stdev','jpg');
    
    figure;
    errorbar(binC{1},betaF_MAP,betaF_MH_std,'-r');
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\betaF(\psi_1) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_bars','fig');
    saveas(gcf,'betaF_eb_bars','jpg');
    
    % - distributions
    figure;

    n_betaFbins = 100;
    betaFbins = linspace(min(betaF_MH(:)),max(betaF_MH(:)),n_betaFbins);
    [X,Y] = meshgrid(binC{1},betaFbins);
    Z = nan(M,n_betaFbins);
    for ii=1:M
        [count,~] = hist(betaF_MH(:,ii),betaFbins);
        Z(ii,:) = count / ( sum(count) * (betaFbins(2)-betaFbins(1)) );
    end

    n_contours = 50;
    contour3(X,Y,Z',n_contours)
    colorbar;
    hold on
    plot3(binC{1},betaF_MAP,zeros(size(betaF_MAP)),'-r');
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\betaF(\psi_1) / -','fontsize',fsize);
    zlabel('pdf(\betaF; \psi_1) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_dist','fig');
    saveas(gcf,'betaF_eb_dist','jpg');
     
elseif dim==2
    
    % pdf
    
    % - naked
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(pdf_MAP,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('pdf(\psi_1,\psi_2) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf','fig');
    saveas(gcf,'pdf','jpg');
    
    % - traces
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(pdf_MAP,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    hold on
    for ii=1:nSamples_MH
        Z = reshape(pdf_MH(ii,:),length(binC{2}),length(binC{1}));
        mesh(X,Y,Z,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
    end
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('pdf(\psi_1,\psi_2) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_traces','fig');
    saveas(gcf,'pdf_traces','jpg');
    
    % - stdev errorbars
    pdf_MH_std = std(pdf_MH,0,1);
    
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(pdf_MAP,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    hold on
    Z = reshape(pdf_MAP+pdf_MH_std,length(binC{2}),length(binC{1}));
    mesh(X,Y,Z,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
    Z = reshape(pdf_MAP-pdf_MH_std,length(binC{2}),length(binC{1}));
    mesh(X,Y,Z,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('pdf(\psi_1,\psi_2) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_eb_limits','fig');
    saveas(gcf,'pdf_eb_limits','jpg');
    
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(pdf_MH_std,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('stdev(pdf(\psi_1,\psi_2)) / (a.u.)^{-1}','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_eb_stdev','fig');
    saveas(gcf,'pdf_eb_stdev','jpg');
    
    
    % betaF
    
    % - naked
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(betaF_MAP,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\betaF(\psi_1,\psi_2) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF','fig');
    saveas(gcf,'betaF','jpg');
    
    % - traces
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(betaF_MAP,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    hold on
    for ii=1:nSamples_MH
        Z = reshape(betaF_MH(ii,:),length(binC{2}),length(binC{1}));
        mesh(X,Y,Z,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
    end
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\betaF(\psi_1,\psi_2) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_traces','fig');
    saveas(gcf,'betaF_traces','jpg');
    
    % - stdev errorbars
    betaF_MH_std = std(betaF_MH,0,1);
    
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(betaF_MAP,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    hold on
    Z = reshape(betaF_MAP+betaF_MH_std,length(binC{2}),length(binC{1}));
    mesh(X,Y,Z,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
    Z = reshape(betaF_MAP-betaF_MH_std,length(binC{2}),length(binC{1}));
    mesh(X,Y,Z,'FaceColor','none','EdgeColor',[0.5 0.5 0.5]);
    hold off
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\betaF(\psi_1,\psi_2) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_limits','fig');
    saveas(gcf,'betaF_eb_limits','jpg');
    
    figure;
    [X,Y] = meshgrid(binC{1},binC{2});
    Z = reshape(betaF_MH_std,length(binC{2}),length(binC{1}));
    meshc(X,Y,Z)
    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('stdev(\betaF(\psi_1,\psi_2)) / -','fontsize',fsize);
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_stdev','fig');
    saveas(gcf,'betaF_eb_stdev','jpg');
    
elseif dim==3

    % pdf
    
    % - naked
    figure;
    
    [X,Y,Z] = meshgrid(binC{1},binC{2},binC{3});
    F = permute(reshape(pdf_MAP,length(binC{3}),length(binC{2}),length(binC{1})),[3,2,1]);
    F = permute(F,[2,1,3]);
    cmin = min(F(~isnan(F)));
    cmax = max(F(~isnan(F)));
    
    %fprintf('pdf:\n');
    %fprintf('How many evenly spaced isosurfaces over the range [%e,%e] would you like to plot? ',cmin,cmax);
    %n_isos = input('');
    n_isos = 7;
    isos = linspace(cmin,cmax,n_isos);
    
    colormap('default');
    cmap = colormap;
    cmap_size = size(cmap);
    cmap_rows = cmap_size(1);

    hold on
    for i=1:n_isos
        
        % interpolating for isosurface color according to colormap
        col_interp = (isos(i) - cmin) / (cmax - cmin) * (cmap_rows-1) + 1;
        col_vec_floor = cmap(floor(col_interp),:);
        col_vec_ceil = cmap(ceil(col_interp),:);
        col_vec = mean(cat(1,col_vec_floor,col_vec_ceil));

        p = patch(isosurface(X,Y,Z,F,isos(i)));
        set(p,'FaceColor',col_vec,'EdgeColor','none','FaceAlpha',0.2);
        
    end;
    hold off 

    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\psi_3 / a.u.','fontsize',fsize);
    set(gca,'clim',[cmin,cmax]);
    colorbar;
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf','fig');
    saveas(gcf,'pdf','jpg');
    
    % - stdev
    pdf_MH_std = std(pdf_MH,0,1);
    
    figure;
    
    [X,Y,Z] = meshgrid(binC{1},binC{2},binC{3});
    F = permute(reshape(pdf_MH_std,length(binC{3}),length(binC{2}),length(binC{1})),[3,2,1]);
    F = permute(F,[2,1,3]);
    cmin = min(F(~isnan(F)));
    cmax = max(F(~isnan(F)));
    
    %fprintf('std(pdf):\n');
    %fprintf('How many evenly spaced isosurfaces over the range [%e,%e] would you like to plot? ',cmin,cmax);
    %n_isos = input('');
    n_isos = 7;
    isos = linspace(cmin,cmax,n_isos);
    
    colormap('default');
    cmap = colormap;
    cmap_size = size(cmap);
    cmap_rows = cmap_size(1);

    hold on
    for i=1:n_isos
        
        % interpolating for isosurface color according to colormap
        col_interp = (isos(i) - cmin) / (cmax - cmin) * (cmap_rows-1) + 1;
        col_vec_floor = cmap(floor(col_interp),:);
        col_vec_ceil = cmap(ceil(col_interp),:);
        col_vec = mean(cat(1,col_vec_floor,col_vec_ceil));

        p = patch(isosurface(X,Y,Z,F,isos(i)));
        set(p,'FaceColor',col_vec,'EdgeColor','none','FaceAlpha',0.2);
        
    end;
    hold off 

    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\psi_3 / a.u.','fontsize',fsize);
    set(gca,'clim',[cmin,cmax]);
    colorbar;
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'pdf_eb_stdev','fig');
    saveas(gcf,'pdf_eb_stdev','jpg');
    
    
    % betaF
    
    % - naked
    figure;
    
    [X,Y,Z] = meshgrid(binC{1},binC{2},binC{3});
    F = permute(reshape(betaF_MAP,length(binC{3}),length(binC{2}),length(binC{1})),[3,2,1]);
    F = permute(F,[2,1,3]);
    cmin = min(F(~isnan(F)));
    cmax = max(F(~isnan(F)));
    
    %fprintf('betaF:\n');
    %fprintf('How many evenly spaced isosurfaces over the range [%e,%e] would you like to plot? ',cmin,cmax);
    %n_isos = input('');
    n_isos = 7;
    isos = linspace(cmin,cmax,n_isos);
    
    colormap('default');
    cmap = colormap;
    cmap_size = size(cmap);
    cmap_rows = cmap_size(1);

    hold on
    for i=1:n_isos
        
        % interpolating for isosurface color according to colormap
        col_interp = (isos(i) - cmin) / (cmax - cmin) * (cmap_rows-1) + 1;
        col_vec_floor = cmap(floor(col_interp),:);
        col_vec_ceil = cmap(ceil(col_interp),:);
        col_vec = mean(cat(1,col_vec_floor,col_vec_ceil));

        p = patch(isosurface(X,Y,Z,F,isos(i)));
        set(p,'FaceColor',col_vec,'EdgeColor','none','FaceAlpha',0.2);
        
    end;
    hold off 

    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\psi_3 / a.u.','fontsize',fsize);
    set(gca,'clim',[cmin,cmax]);
    colorbar;
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF','fig');
    saveas(gcf,'betaF','jpg');
    
    % - stdev
    betaF_MH_std = std(betaF_MH,0,1);
    
    figure;
    
    [X,Y,Z] = meshgrid(binC{1},binC{2},binC{3});
    F = permute(reshape(betaF_MH_std,length(binC{3}),length(binC{2}),length(binC{1})),[3,2,1]);
    F = permute(F,[2,1,3]);
    cmin = min(F(~isnan(F)));
    cmax = max(F(~isnan(F)));
    
    %fprintf('std(betaF):\n');
    %fprintf('How many evenly spaced isosurfaces over the range [%e,%e] would you like to plot? ',cmin,cmax);
    %n_isos = input('');
    n_isos = 7;
    isos = linspace(cmin,cmax,n_isos);
    
    colormap('default');
    cmap = colormap;
    cmap_size = size(cmap);
    cmap_rows = cmap_size(1);

    hold on
    for i=1:n_isos
        
        % interpolating for isosurface color according to colormap
        col_interp = (isos(i) - cmin) / (cmax - cmin) * (cmap_rows-1) + 1;
        col_vec_floor = cmap(floor(col_interp),:);
        col_vec_ceil = cmap(ceil(col_interp),:);
        col_vec = mean(cat(1,col_vec_floor,col_vec_ceil));

        p = patch(isosurface(X,Y,Z,F,isos(i)));
        set(p,'FaceColor',col_vec,'EdgeColor','none','FaceAlpha',0.2);
        
    end;
    hold off 

    xlabel('\psi_1 / a.u.','fontsize',fsize);
    ylabel('\psi_2 / a.u.','fontsize',fsize);
    zlabel('\psi_3 / a.u.','fontsize',fsize);
    set(gca,'clim',[cmin,cmax]);
    colorbar;
    set(gca,'fontsize',fsize);
    set(gcf,'color','w');
    saveas(gcf,'betaF_eb_stdev','fig');
    saveas(gcf,'betaF_eb_stdev','jpg');
    
else
    
    error('Dimensionality of data dim = %d; plotting only available for {1,2,3}-dimensional space',dim)
    
end

fprintf('DONE!\n\n')



end
