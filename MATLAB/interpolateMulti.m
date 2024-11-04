function OUT = interpolateMulti(IN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code for multiscale kernel interpolation
%%% Bermanis et al. Appl Comput Harm Anal 2013 [exact matching version, applied for regression]
%%% Duchateau et al. SEE Geom Sci Inf 2013 [inexact matching version, applied to map data to a manifold]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Looking for y = f(x) with known samples yi = f(xi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading data and options
data.x = IN.xi;
data.y = IN.yi;
data.numS = size(data.x,1);
data.dim  = size(data.y,2);
testData.x = IN.x;
testData.numS = size(testData.x,1);

if isfield(IN,'kNN') %%% number of nearest neighbors used to estimate the density of samples
    kNN = IN.kNN;
else
    kNN = 10;
end
if isfield(IN,'singleScale') %%% use a single scale
    useSingleScale = 1;
    data.T = IN.singleScale;
else
    useSingleScale = 0;
    data.T = NaN;
end
if isfield(IN,'densityFactor') %%% at which scale to stop the iterations (density * factor)
    densityFactor = IN.densityFactor;
else
    densityFactor = 1;
end
if isfield(IN,'gammaList') %%% weight(s) between regularization and similarity (the lower the smoother)
    data.gammaList = IN.gammaList;
else
    data.gammaList = 10.^0;
end
if isfield(IN,'usePINV') %%% use pinv for inverse computations
    usePINV = IN.usePINV;
else
    usePINV = 0;
end
if isfield(IN,'useNYSTROM') %%% use the definition of Nystrom extension
    useNYSTROM = IN.useNYSTROM;
else
    useNYSTROM = 0;
end
if isfield(IN,'exactPoints')
    exactPoints = IN.exactPoints; %%% list of samples where to have exact matching
else
    exactPoints = [];
end

clear IN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dg = squareform(pdist(data.x));

if(useSingleScale == 0)
    %% Estimating data density
    diam = max(Dg(:));
    data.T = diam^2;
    tmp = Dg + diag( Inf(data.numS,1) );  %% put diagonal coefficients to -1
    tmpB = sort(tmp,1,'ascend');
    clear tmp;
    tmpB = tmpB(1:min(kNN,data.numS-1),:); %% approximate density from k neighbors
    tmpB = mean(tmpB,1);
    data.density = mean(tmpB(:));
    clear tmpB;

    it_max = 20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiscale regression - inexact (terminology was kept similar to Bermanis et al. Appl Comput Harm Anal 2013)

x = data.x;
x_star = testData.x;
f = data.y;

OUT.OUT = cell(length(data.gammaList),1);
for gI = 1:length(data.gammaList)
    disp(['gI #',num2str(gI)]);

    gammaI = data.gammaList(gI);

    gamma = gammaI;
    F_s_old = zeros(data.numS,data.dim);
    F_star_s_old = zeros(testData.numS,data.dim);
    
    s = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(useSingleScale == 1)
        e_s = data.T / 2^s;
        Ge_s = exp( -Dg.^2 / e_s );  %% (n,n)
        pointsInexact = 1:size(Ge_s,1);
        pointsInexact(exactPoints) = [];
        output_03 = MINEXACT_Bermanis_ACHA_2013_ALGO_03( Ge_s , f - F_s_old , x , x_star , e_s , gamma , pointsInexact , useNYSTROM , usePINV);
        F_star_s_old = output_03.f_star_s;
    else
        %%% Stop when the kernel size is smaller than the density of samples
        while ( ( s <= it_max ) && ( sqrt(data.T/2^s) > (densityFactor*data.density) ) )
            fprintf([num2str(s),' ']);
            e_s = data.T / 2^s;
            Ge_s = exp( -Dg.^2 / e_s );  %% (n,n)
            pointsInexact = 1:size(Ge_s,1);
            pointsInexact(exactPoints) = [];
            output_03 = MINEXACT_Bermanis_ACHA_2013_ALGO_03( Ge_s , f - F_s_old , x , x_star , e_s , gamma , pointsInexact , useNYSTROM , usePINV);
            F_s = F_s_old + output_03.f_s;
            F_star_s = F_star_s_old + output_03.f_star_s;
            s = s+1;
            F_s_old = F_s;
            F_star_s_old = F_star_s;
        end
        OUT.density = data.density;
        OUT.finalT{gI} = [sqrt(data.T/2^(s-1)),s-1];
    end
    
    OUT.OUT{gI} = F_star_s_old;
    fprintf('\n');
end

OUT.diam = sqrt(data.T);
