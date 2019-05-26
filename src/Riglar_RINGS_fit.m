function fitInfo = Riglar_RINGS_fit(im, varargin)

% this function fits the phase profile from Riglar et al,  ringed bacteria
% data in 2 colors (YFP and CFP output channels)


%% parse inputs

% lengths are in millimeters

ip = inputParser;
ip.addParameter('imgCal', 34.7, @(x)(numel(x)==1 && x > 0 )); % calibration of the image: px / mm
ip.addParameter('rMin', 0, @(x)(numel(x)==1 && x >= 0 )); % mask the colony at radius less than this.  in mm's.
ip.addParameter('rMax', 1.3, @(x)(numel(x)==1 && x > 0 )); % mask the colony at radius larger than this.  in mm's. For LPT320 & PAS715 use 1.3, for PAS716 use 1.1
ip.addParameter('expectedSlope', 0.39, @(x)(numel(x)==1 && x > 0 )); % expected slope of the phase profile.  This is in [rad]/[px]. For LPT320 -0.39; PAS715- 0.34; PAS716 - 0.43 
% You can estimate expectedSlope by measuring distance between two consecutive peaks in the pattern.
ip.addParameter('colorPhaseShift', 1.5); % CFP lags YFP by this amount.  In radians. For LPT320 - 1.5; PAS715 - 0.9, PAS716 - 1.0.
ip.addParameter('slopeTol', 0.3, @(x)(numel(x)==1 && x > 0 )); % fractional change in slope of phase profile that is tolerable (outside of this range, the sample will be discarded) For LPT320 0.3; PAS715/716 - 0.4.
ip.addParameter('uncertaintyTol', 0.1, @(x)(numel(x)==1 && x > 0 )); % tolerable uncertainty in the fit of the slope relative to the value of the slope (outside of this range, the sample will be discarded)
ip.addParameter('centerTol', 0.03, @(x)(numel(x)==1 && x > 0 )); % distance that center of colony can be off from center of image, as a fraction of the image width (outside of this range, the sample will be discarded)
ip.addParameter('modelOrder', 1); % set the order of the polynomial describing the phase profile
ip.addParameter('decayModelOrder', 2, @(x)(x==2 || x==4)); % set the order of the polynomial describing the intensity decay profile.  must be 2 or 4.
ip.addParameter('plotFlag', 0); % Set to 1 to output images of model fit.  careful, this may be a lot of figures.

ip.addParameter('residualWeightingMethod', 0);
% 0: no weighting (as before)
% 1: spatial weighting
% 2: robust weighting (down weight outliers compared to standard squared error)
ip.addParameter('rMed', 0.3, @(x)(numel(x)==1 && x > 0 )); % if you select spatial weighting, the weights will be constant within this inner region

ip.parse(varargin{:});

in_params = ip.Results;

%% 

% convert im to double for fitting
im = double(im);

% dims
width  = size(im,2);
height = size(im,1);

%% map of radius and spatial weights

% map of radius
[X, Y] = meshgrid(1:width,1:height);

% duplicate for two colors
X = repmat(X,[1,1,2]);
Y = repmat(Y,[1,1,2]);

center = [(width+1)/2, (height+1)/2];
r = sqrt((X-center(1)).^2 + (Y-center(2)).^2);

% generate corresponding weights
rMedPx = in_params.rMed * in_params.imgCal;
fitWeights = spatialWeightingFnc(r, rMedPx); % mask later, when you have inner and outer radii masks

%% normalize radial intensity profile (not the rings)

% generate mask based on max radius
rMaxPx = in_params.rMax * in_params.imgCal;
mask = r <= rMaxPx;
num_px_per_image = sum(sum(mask(:,:,1))); % for splitting images apart

% apply mask to data for fitting radial intensity decay
imData = im(mask);
xData = [X(mask), Y(mask)];

% fit radially decaying intensity profile
f = @(p) intDecayProfile_2color(p,xData,in_params.decayModelOrder) - imData;

if in_params.decayModelOrder == 2
    pInit = [center(1), center(2), ...
        max(imData(1:num_px_per_image)), max(imData(num_px_per_image+1:end)), ...
        -0.01, -0.01];
    pLB = [1, 1, ...
        min(imData(1:num_px_per_image)), min(imData(num_px_per_image+1:end)), ...
        -1/eps, -1/eps];
    pUB = [width, height, ...
        max(imData(1:num_px_per_image)), max(imData(num_px_per_image+1:end)), ...
        0, 0];
elseif in_params.decayModelOrder == 4
    pInit = [center(1), center(2), ...
        max(imData(1:num_px_per_image)), max(imData(num_px_per_image+1:end)), ...
        -0.01, -0.01, ...
        0, 0];
    pLB = [1, 1, ...
        min(imData(1:num_px_per_image)), min(imData(num_px_per_image+1:end)), ...
        -1/eps, -1/eps, ...
        -1/eps, -1/eps];
    pUB = [width, height, ...
        max(imData(1:num_px_per_image)), max(imData(num_px_per_image+1:end)), ...
        0, 0, ...
        0, 0];
end

pFit = lsqnonlin(f,pInit,pLB,pUB);

% evaluate on whole image
yPred = intDecayProfile_2color(pFit,[X(:),Y(:)],in_params.decayModelOrder);
yPred = reshape(yPred, [height,width,2]);

% normalize data for fitting rings
normalized = reshape(im ./ yPred, [height,width,2]);

% apply mask
yPred(mask==0) = NaN;
normalized(mask==0) = NaN;
imMasked = im;
imMasked(mask==0) = NaN;

% to visualize radial profile, uncomment this section
% if in_params.plotFlag
%     
%     figure
%     % yfp
%     subplot(2,3,1)
%     imagesc(imMasked(:,:,1))
%     subplot(2,3,2)
%     imagesc(yPred(:,:,1))
%     subplot(2,3,3)
%     imagesc(normalized(:,:,1))
% 
%     % cfp
%     subplot(2,3,4)
%     imagesc(imMasked(:,:,2))
%     subplot(2,3,5)
%     imagesc(yPred(:,:,2))
%     subplot(2,3,6)
%     imagesc(normalized(:,:,2))

% end

%% set up optimization problem for fitting rings

% user params
maxIter = 5e2;
% tolFun
% tolX

% set up for different types of weights
switch in_params.residualWeightingMethod
    case 0
        opts = statset('MaxIter', maxIter);
    case 1
        opts = statset('MaxIter', maxIter, 'RobustWgtFun', []);
    case 2
        opts = statset('MaxIter', maxIter, 'RobustWgtFun', 'cauchy');
    otherwise
        error('fit weights not supported')
end

%% fit normalized with radial sinusoid

% generate mask based on min and max radii
rMinPx = in_params.rMin * in_params.imgCal;
rMaxPx = in_params.rMax * in_params.imgCal;
innerMask = r > rMinPx;
outerMask = r <= rMaxPx;
mask = logical(innerMask .* outerMask);
num_px_per_image = sum(sum(mask(:,:,1))); % for splitting images apart

% apply new mask to data for fitting rings
imData = normalized(mask);
xData = [X(mask), Y(mask)];

% scale parameters to avoid optimization problems
pScaling = [center(1), center(2), ...
    max(imData(1:num_px_per_image)), 0.1, ...
    max(imData(num_px_per_image+1:end)), 0.1, ...
    pi, in_params.expectedSlope, ones(1,in_params.modelOrder-1)] .^ -1;

% fit function
fitFun = @(q,xData) phaseProfile_2color(q,xData,in_params.modelOrder,pScaling,in_params.colorPhaseShift);

% init params
qInit = [center(1), center(2), ...
    max(imData(1:num_px_per_image)), 0.1, ...
    max(imData(num_px_per_image+1:end)), 0.1, ...
    0, in_params.expectedSlope, zeros(1,in_params.modelOrder-1)];

% fit
switch in_params.residualWeightingMethod
    case 1
        fitWeights = fitWeights(mask);
        [qFit, R, J, CovQ, ~] = nlinfit(xData, imData, fitFun, qInit.* pScaling, opts, 'Weights', fitWeights);
    otherwise
        [qFit, R, J, CovQ, ~] = nlinfit(xData, imData, fitFun, qInit.* pScaling, opts);
end

fitInfo.Res = R;
fitInfo.J = full(J);
fitInfo.Cov = CovQ;

%% Parameter covariance and CI estimates

% Undo this scaling last to avoid numerical issues
qFit = qFit ./ pScaling;
fitInfo.J = bsxfun(@rdivide,fitInfo.J,pScaling);
fitInfo.Cov = fitInfo.Cov ./ (pScaling' * pScaling); % Undo scaling
fitInfo.ParScaling = pScaling;

% Parameter confidence intervals
pCi = nlparci(qFit, fitInfo.Res, 'covar', fitInfo.Cov);

% Put it all in a structure for output
extra_in_params = {'second_order_term','third_order_term','fourth_order_term','fifth_order_term'};
fitInfo.ParOrder = [{'center_X','center_Y','intensity_offset_yfp','intensity_amplitude_yfp','intensity_offset_cfp','intensity_amplitude_cfp','theta_0','slope'} extra_in_params(1:in_params.modelOrder-1)];

for j = 1:numel(qFit)
    fitInfo.(fitInfo.ParOrder{j}) = qFit(j);
    fitInfo.([fitInfo.ParOrder{j} '_CI']) = pCi(j,:);
end

%% correct for negative amplitude (if it ever happens)

if (fitInfo.intensity_amplitude_yfp < 0) && (fitInfo.intensity_amplitude_cfp < 0)
    fitInfo.intensity_amplitude_yfp = abs(fitInfo.intensity_amplitude_yfp);
    fitInfo.intensity_amplitude_cfp = abs(fitInfo.intensity_amplitude_cfp);
    fitInfo.theta_0 = fitInfo.theta_0 + pi;
end

%% Quality of fit metric

% do parameters lie outside acceptable range?
center_offset = sqrt((qFit(1)-center(1))^2 + (qFit(2)-center(2))^2);
slope_offset  = abs(qFit(8) - in_params.expectedSlope);

% estimate uncertainty in slope
slope_unc = sqrt(fitInfo.Cov(6,6));

fitInfo.qualityFlag = 0;

% set quality flag
if center_offset/width > in_params.centerTol
    warning('center of colony too far from center of image')
    fitInfo.qualityFlag = 2;
elseif slope_offset/in_params.expectedSlope > in_params.slopeTol
    warning('slope is not within tolerance')
    fitInfo.qualityFlag = 3;    
elseif slope_unc/qFit(8) > in_params.uncertaintyTol
    warning('uncertainty in fitted slope is not within tolerance')
    fitInfo.qualityFlag = 4;
elseif sign(fitInfo.intensity_amplitude_yfp) ~= sign(fitInfo.intensity_amplitude_cfp)
    warning('cfp and yfp amplitude have opposite sign.  fit is unrealiable')
    fitInfo.qualityFlag = 5;
end

%% visualize fit

if in_params.plotFlag
    
    % make images from predictions
    yPred = phaseProfile_2color(qFit,[X(:), Y(:)],in_params.modelOrder,ones(size(qFit)),in_params.colorPhaseShift);
    yPred = reshape(yPred, [height, width, 2]);
    
    yRes = yPred - normalized;

    % apply mask
    yPred(mask==0) = NaN;
    yRes(mask==0) = NaN;
    imMasked = normalized;
    imMasked(mask==0) = NaN;
    
    figure
    subplot(2,3,1)
    imagesc(imMasked(:,:,1))
    subplot(2,3,2)
    imagesc(yPred(:,:,1))
    title(['Quality Flag = ', num2str(fitInfo.qualityFlag)])
    subplot(2,3,3)
    imagesc(yRes(:,:,1))
    
    subplot(2,3,4)
    imagesc(imMasked(:,:,2))
    subplot(2,3,5)
    imagesc(yPred(:,:,2))
    title(['Quality Flag = ', num2str(fitInfo.qualityFlag)])
    subplot(2,3,6)
    imagesc(yRes(:,:,2))

end

end

%% functions for fitting intensity profile

function y = intDecayProfile_2color(p,xData,modelOrder)

    xc = p(1);
    yc = p(2);
    r = sqrt((xData(:,1) - xc).^2 + (xData(:,2) - yc).^2);
    
    % split data from two separate colors
    assert(mod(numel(r),2)==0)
    num_px_per_img = numel(r)/2;
    r1 = r(1:num_px_per_img);
    r2 = r(num_px_per_img+1:end);
    assert(isequal(r1,r2))
    r = r1;

    switch modelOrder
        case 2
            y1 = p(3) + p(5)*r.^2;
            y2 = p(4) + p(6)*r.^2;
        case 4
            y1 = p(3) + p(5)*r.^2 + p(7)*r.^4;
            y2 = p(4) + p(6)*r.^2 + p(8)*r.^4;
        otherwise
            error('decay model order')
    end
    
    y = [y1;y2];
    
end

function y = phaseProfile_2color(q,xData,modelOrder,pScaling,phaseShift)

    % scale parameters
    q = q ./ pScaling;
    
    % radius for every pixel
    xc = q(1);
    yc = q(2);
    r = sqrt((xData(:,1) - xc).^2 + (xData(:,2) - yc).^2);
    
    % split data from two separate colors
    assert(mod(numel(r),2)==0)
    num_px_per_img = numel(r)/2;
    r1 = r(1:num_px_per_img);
    r2 = r(num_px_per_img+1:end);
    assert(isequal(r1,r2))
    r = r1;

    % phase
    switch modelOrder
        case 1
            phase1 = phaseShift + q(7) + q(8)*r;
            phase2 = q(7) + q(8)*r;
        case 2
            phase1 = phaseShift + q(7) + q(8)*r + q(9)*(r.^2);
            phase2 = q(7) + q(8)*r + q(9)*(r.^2);
        case 3
            phase1 = phaseShift + q(7) + q(8)*r + q(9)*(r.^2) + q(10)*(r.^3);
            phase2 = q(7) + q(8)*r + q(9)*(r.^2) + q(10)*(r.^3);
        otherwise
            error('model order too high')
    end
    
    % model output
    y1 = q(3) + q(4)*sin(phase1);
    y2 = q(5) + q(6)*sin(phase2);
    
    y = [y1;y2];
    
end

function weights = spatialWeightingFnc(r,rMed)

    weights = (1/rMed)*ones(size(r));
    
    % outer region with decaying weights
    mask = r > rMed;
    weights(mask) = 1./r(mask);
    
end