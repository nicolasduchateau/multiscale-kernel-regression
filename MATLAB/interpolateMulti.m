{"payload":{"allShortcutsEnabled":true,"fileTree":{"":{"items":[{"name":"MATLAB","path":"MATLAB","contentType":"directory"},{"name":"LICENSE.txt","path":"LICENSE.txt","contentType":"file"},{"name":"MINEXACT_Bermanis_ACHA_2013_ALGO_03.m","path":"MINEXACT_Bermanis_ACHA_2013_ALGO_03.m","contentType":"file"},{"name":"README.md","path":"README.md","contentType":"file"},{"name":"interpolateMulti.m","path":"interpolateMulti.m","contentType":"file"},{"name":"testSIN.m","path":"testSIN.m","contentType":"file"}],"totalCount":6}},"fileTreeProcessingTime":1.809169,"foldersToFetch":[],"reducedMotionEnabled":"system","repo":{"id":529246634,"defaultBranch":"main","name":"multiscale-kernel-regression","ownerLogin":"nicolasduchateau","currentUserCanPush":true,"isFork":false,"isEmpty":false,"createdAt":"2022-08-26T14:21:21.000+02:00","ownerAvatar":"https://avatars.githubusercontent.com/u/30652123?v=4","public":true,"private":false,"isOrgOwned":false},"symbolsExpanded":false,"treeExpanded":true,"refInfo":{"name":"main","listCacheKey":"v0:1663866967.367121","canEdit":true,"refType":"branch","currentOid":"6011fbb9f792c873099d200850490cf05200a55a"},"path":"interpolateMulti.m","currentUser":{"id":30652123,"login":"nicolasduchateau","userEmail":"n.duchateaucistib@gmail.com"},"blob":{"rawLines":["function OUT = interpolateMulti(IN)","","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%%% Code for multiscale kernel interpolation","%%% Bermanis et al. Appl Comput Harm Anal 2013 [exact matching version, applied for regression]","%%% Duchateau et al. SEE Geom Sci Inf 2013 [inexact matching version, applied to map data to a manifold]","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","","%%% Looking for y = f(x) with known samples yi = f(xi)","","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%% Loading data and options","data.x = IN.xi;","data.y = IN.yi;","data.numS = size(data.x,1);","data.dim  = size(data.y,2);","testData.x = IN.x;","testData.numS = size(testData.x,1);","","if isfield(IN,'kNN') %%% number of nearest neighbors used to estimate the density of samples","    kNN = IN.kNN;","else","    kNN = 10;","end","if isfield(IN,'singleScale') %%% use a single scale","    useSingleScale = 1;","    data.T = IN.singleScale;","else","    useSingleScale = 0;","    data.T = NaN;","end","if isfield(IN,'densityFactor') %%% at which scale to stop the iterations (density * factor)","    densityFactor = IN.densityFactor;","else","    densityFactor = 1;","end","if isfield(IN,'gammaList') %%% weight(s) between regularization and similarity (the lower the smoother)","    data.gammaList = IN.gammaList;","else","    data.gammaList = 10.^0;","end","if isfield(IN,'usePINV') %%% use pinv for inverse computations","    usePINV = IN.usePINV;","else","    usePINV = 0;","end","if isfield(IN,'useNYSTROM') %%% use the definition of Nystrom extension","    useNYSTROM = IN.useNYSTROM;","else","    useNYSTROM = 0;","end","if isfield(IN,'exactPoints')","    exactPoints = IN.exactPoints; %%% list of samples where to have exact matching","else","    exactPoints = [];","end","","clear IN;","","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","Dg = squareform(pdist(data.x));","","if(useSingleScale == 0)","    %% Estimating data density","    diam = max(Dg(:));","    data.T = diam^2;","    tmp = Dg + diag( Inf(data.numS,1) );  %% put diagonal coefficients to -1","    tmpB = sort(tmp,1,'ascend');","    clear tmp;","    tmpB = tmpB(1:min(kNN,data.numS-1),:); %% approximate density from k neighbors","    tmpB = mean(tmpB,1);","    data.density = mean(tmpB(:));","    clear tmpB;","","    it_max = 20;","end","","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","%% Multiscale regression - inexact (terminology was kept similar to Bermanis et al. Appl Comput Harm Anal 2013)","","x = data.x;","x_star = testData.x;","f = data.y;","","OUT.OUT = cell(length(data.gammaList),1);","for gI = 1:length(data.gammaList)","    disp(['gI #',num2str(gI)]);","","    gammaI = data.gammaList(gI);","","    gamma = gammaI;","    F_s_old = zeros(data.numS,data.dim);","    F_star_s_old = zeros(testData.numS,data.dim);","    ","    s = 0;","    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%","    if(useSingleScale == 1)","        e_s = data.T / 2^s;","        Ge_s = exp( -Dg.^2 / e_s );  %% (n,n)","        pointsInexact = 1:size(Ge_s,1);","        pointsInexact(exactPoints) = [];","        output_03 = MINEXACT_Bermanis_ACHA_2013_ALGO_03( Ge_s , f - F_s_old , x , x_star , e_s , gamma , pointsInexact , useNYSTROM , usePINV);","        F_star_s_old = output_03.f_star_s;","    else","        %%% Stop when the kernel size is smaller than the density of samples","        while ( ( s <= it_max ) && ( sqrt(data.T/2^s) > (densityFactor*data.density) ) )","            fprintf([num2str(s),' ']);","            e_s = data.T / 2^s;","            Ge_s = exp( -Dg.^2 / e_s );  %% (n,n)","            pointsInexact = 1:size(Ge_s,1);","            pointsInexact(exactPoints) = [];","            output_03 = MINEXACT_Bermanis_ACHA_2013_ALGO_03( Ge_s , f - F_s_old , x , x_star , e_s , gamma , pointsInexact , useNYSTROM , usePINV);","            F_s = F_s_old + output_03.f_s;","            F_star_s = F_star_s_old + output_03.f_star_s;","            s = s+1;","            F_s_old = F_s;","            F_star_s_old = F_star_s;","        end","        OUT.density = data.density;","        OUT.finalT{gI} = [sqrt(data.T/2^(s-1)),s-1];","    end","    ","    OUT.OUT{gI} = F_star_s_old;","    fprintf('\\n');","end","","OUT.diam = sqrt(data.T);"],"stylingDirectives":[[{"start":0,"end":8,"cssClass":"pl-k"},{"start":9,"end":12,"cssClass":"pl-v"},{"start":13,"end":14,"cssClass":"pl-k"},{"start":15,"end":31,"cssClass":"pl-en"},{"start":32,"end":34,"cssClass":"pl-v"}],[],[{"start":0,"end":98,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":98,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":44,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":95,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":104,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":98,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":98,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[],[{"start":0,"end":54,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[],[{"start":0,"end":98,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":27,"cssClass":"pl-c"},{"start":0,"end":2,"cssClass":"pl-c"},{"start":3,"end":27,"cssClass":"pl-en"}],[{"start":7,"end":8,"cssClass":"pl-k"},{"start":9,"end":11,"cssClass":"pl-smi"},{"start":12,"end":14,"cssClass":"pl-smi"}],[{"start":7,"end":8,"cssClass":"pl-k"},{"start":9,"end":11,"cssClass":"pl-smi"},{"start":12,"end":14,"cssClass":"pl-smi"}],[{"start":10,"end":11,"cssClass":"pl-k"},{"start":12,"end":16,"cssClass":"pl-en"},{"start":17,"end":21,"cssClass":"pl-smi"},{"start":22,"end":23,"cssClass":"pl-smi"},{"start":24,"end":25,"cssClass":"pl-c1"}],[{"start":10,"end":11,"cssClass":"pl-k"},{"start":12,"end":16,"cssClass":"pl-en"},{"start":17,"end":21,"cssClass":"pl-smi"},{"start":22,"end":23,"cssClass":"pl-smi"},{"start":24,"end":25,"cssClass":"pl-c1"}],[{"start":11,"end":12,"cssClass":"pl-k"},{"start":13,"end":15,"cssClass":"pl-smi"},{"start":16,"end":17,"cssClass":"pl-smi"}],[{"start":14,"end":15,"cssClass":"pl-k"},{"start":16,"end":20,"cssClass":"pl-en"},{"start":21,"end":29,"cssClass":"pl-smi"},{"start":30,"end":31,"cssClass":"pl-smi"},{"start":32,"end":33,"cssClass":"pl-c1"}],[],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":19,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":18,"end":19,"cssClass":"pl-pds"},{"start":21,"end":92,"cssClass":"pl-c"},{"start":21,"end":22,"cssClass":"pl-c"}],[{"start":8,"end":9,"cssClass":"pl-k"},{"start":10,"end":12,"cssClass":"pl-smi"},{"start":13,"end":16,"cssClass":"pl-smi"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":8,"end":9,"cssClass":"pl-k"},{"start":10,"end":12,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":27,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":26,"end":27,"cssClass":"pl-pds"},{"start":29,"end":51,"cssClass":"pl-c"},{"start":29,"end":30,"cssClass":"pl-c"}],[{"start":19,"end":20,"cssClass":"pl-k"},{"start":21,"end":22,"cssClass":"pl-c1"}],[{"start":11,"end":12,"cssClass":"pl-k"},{"start":13,"end":15,"cssClass":"pl-smi"},{"start":16,"end":27,"cssClass":"pl-smi"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":19,"end":20,"cssClass":"pl-k"},{"start":21,"end":22,"cssClass":"pl-c1"}],[{"start":11,"end":12,"cssClass":"pl-k"},{"start":13,"end":16,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":29,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":28,"end":29,"cssClass":"pl-pds"},{"start":31,"end":91,"cssClass":"pl-c"},{"start":31,"end":32,"cssClass":"pl-c"}],[{"start":18,"end":19,"cssClass":"pl-k"},{"start":20,"end":22,"cssClass":"pl-smi"},{"start":23,"end":36,"cssClass":"pl-smi"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":18,"end":19,"cssClass":"pl-k"},{"start":20,"end":21,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":25,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":24,"end":25,"cssClass":"pl-pds"},{"start":27,"end":103,"cssClass":"pl-c"},{"start":27,"end":28,"cssClass":"pl-c"}],[{"start":19,"end":20,"cssClass":"pl-k"},{"start":21,"end":23,"cssClass":"pl-smi"},{"start":24,"end":33,"cssClass":"pl-smi"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":19,"end":20,"cssClass":"pl-k"},{"start":21,"end":23,"cssClass":"pl-c1"},{"start":23,"end":25,"cssClass":"pl-k"},{"start":25,"end":26,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":23,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":22,"end":23,"cssClass":"pl-pds"},{"start":25,"end":62,"cssClass":"pl-c"},{"start":25,"end":26,"cssClass":"pl-c"}],[{"start":12,"end":13,"cssClass":"pl-k"},{"start":14,"end":16,"cssClass":"pl-smi"},{"start":17,"end":24,"cssClass":"pl-smi"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":12,"end":13,"cssClass":"pl-k"},{"start":14,"end":15,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":26,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":25,"end":26,"cssClass":"pl-pds"},{"start":28,"end":71,"cssClass":"pl-c"},{"start":28,"end":29,"cssClass":"pl-c"}],[{"start":15,"end":16,"cssClass":"pl-k"},{"start":17,"end":19,"cssClass":"pl-smi"},{"start":20,"end":30,"cssClass":"pl-smi"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":15,"end":16,"cssClass":"pl-k"},{"start":17,"end":18,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":10,"cssClass":"pl-en"},{"start":11,"end":13,"cssClass":"pl-smi"},{"start":14,"end":27,"cssClass":"pl-s"},{"start":14,"end":15,"cssClass":"pl-pds"},{"start":26,"end":27,"cssClass":"pl-pds"}],[{"start":16,"end":17,"cssClass":"pl-k"},{"start":18,"end":20,"cssClass":"pl-smi"},{"start":21,"end":32,"cssClass":"pl-smi"},{"start":34,"end":82,"cssClass":"pl-c"},{"start":34,"end":35,"cssClass":"pl-c"}],[{"start":0,"end":4,"cssClass":"pl-k"}],[{"start":16,"end":17,"cssClass":"pl-k"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[],[{"start":0,"end":5,"cssClass":"pl-en"},{"start":6,"end":8,"cssClass":"pl-s"}],[],[{"start":0,"end":49,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":3,"end":4,"cssClass":"pl-k"},{"start":5,"end":15,"cssClass":"pl-en"},{"start":16,"end":21,"cssClass":"pl-en"},{"start":22,"end":26,"cssClass":"pl-smi"},{"start":27,"end":28,"cssClass":"pl-smi"}],[],[{"start":0,"end":2,"cssClass":"pl-k"},{"start":3,"end":17,"cssClass":"pl-smi"},{"start":18,"end":20,"cssClass":"pl-k"},{"start":21,"end":22,"cssClass":"pl-c1"}],[{"start":4,"end":30,"cssClass":"pl-c"},{"start":4,"end":6,"cssClass":"pl-c"},{"start":7,"end":30,"cssClass":"pl-en"}],[{"start":9,"end":10,"cssClass":"pl-k"},{"start":11,"end":14,"cssClass":"pl-en"},{"start":15,"end":17,"cssClass":"pl-en"},{"start":18,"end":19,"cssClass":"pl-k"}],[{"start":11,"end":12,"cssClass":"pl-k"},{"start":13,"end":17,"cssClass":"pl-smi"},{"start":17,"end":18,"cssClass":"pl-k"},{"start":18,"end":19,"cssClass":"pl-c1"}],[{"start":8,"end":9,"cssClass":"pl-k"},{"start":10,"end":12,"cssClass":"pl-smi"},{"start":13,"end":14,"cssClass":"pl-k"},{"start":15,"end":19,"cssClass":"pl-en"},{"start":21,"end":24,"cssClass":"pl-en"},{"start":21,"end":24,"cssClass":"pl-c1"},{"start":25,"end":29,"cssClass":"pl-smi"},{"start":30,"end":34,"cssClass":"pl-smi"},{"start":35,"end":36,"cssClass":"pl-c1"},{"start":42,"end":76,"cssClass":"pl-c"},{"start":42,"end":44,"cssClass":"pl-c"},{"start":45,"end":76,"cssClass":"pl-en"}],[{"start":9,"end":10,"cssClass":"pl-k"},{"start":11,"end":15,"cssClass":"pl-en"},{"start":16,"end":19,"cssClass":"pl-smi"},{"start":20,"end":21,"cssClass":"pl-c1"},{"start":22,"end":30,"cssClass":"pl-s"},{"start":22,"end":23,"cssClass":"pl-pds"},{"start":29,"end":30,"cssClass":"pl-pds"}],[{"start":4,"end":9,"cssClass":"pl-en"},{"start":10,"end":13,"cssClass":"pl-s"}],[{"start":9,"end":10,"cssClass":"pl-k"},{"start":11,"end":15,"cssClass":"pl-en"},{"start":16,"end":17,"cssClass":"pl-c1"},{"start":17,"end":18,"cssClass":"pl-k"},{"start":18,"end":21,"cssClass":"pl-en"},{"start":22,"end":25,"cssClass":"pl-smi"},{"start":26,"end":30,"cssClass":"pl-smi"},{"start":31,"end":35,"cssClass":"pl-smi"},{"start":35,"end":36,"cssClass":"pl-k"},{"start":36,"end":37,"cssClass":"pl-c1"},{"start":39,"end":40,"cssClass":"pl-k"},{"start":43,"end":82,"cssClass":"pl-c"},{"start":43,"end":45,"cssClass":"pl-c"},{"start":46,"end":82,"cssClass":"pl-en"}],[{"start":9,"end":10,"cssClass":"pl-k"},{"start":11,"end":15,"cssClass":"pl-en"},{"start":16,"end":20,"cssClass":"pl-smi"},{"start":21,"end":22,"cssClass":"pl-c1"}],[{"start":17,"end":18,"cssClass":"pl-k"},{"start":19,"end":23,"cssClass":"pl-en"},{"start":24,"end":28,"cssClass":"pl-en"},{"start":29,"end":30,"cssClass":"pl-k"}],[{"start":4,"end":9,"cssClass":"pl-en"},{"start":10,"end":14,"cssClass":"pl-s"}],[],[{"start":11,"end":12,"cssClass":"pl-k"},{"start":13,"end":15,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[],[{"start":0,"end":49,"cssClass":"pl-c"},{"start":0,"end":1,"cssClass":"pl-c"}],[{"start":0,"end":111,"cssClass":"pl-c"},{"start":0,"end":2,"cssClass":"pl-c"},{"start":3,"end":111,"cssClass":"pl-en"}],[],[{"start":2,"end":3,"cssClass":"pl-k"},{"start":4,"end":8,"cssClass":"pl-smi"},{"start":9,"end":10,"cssClass":"pl-smi"}],[{"start":7,"end":8,"cssClass":"pl-k"},{"start":9,"end":17,"cssClass":"pl-smi"},{"start":18,"end":19,"cssClass":"pl-smi"}],[{"start":2,"end":3,"cssClass":"pl-k"},{"start":4,"end":8,"cssClass":"pl-smi"},{"start":9,"end":10,"cssClass":"pl-smi"}],[],[{"start":8,"end":9,"cssClass":"pl-k"},{"start":10,"end":14,"cssClass":"pl-en"},{"start":15,"end":21,"cssClass":"pl-en"},{"start":22,"end":26,"cssClass":"pl-smi"},{"start":27,"end":36,"cssClass":"pl-smi"},{"start":38,"end":39,"cssClass":"pl-c1"}],[{"start":0,"end":3,"cssClass":"pl-k"},{"start":7,"end":8,"cssClass":"pl-k"},{"start":9,"end":10,"cssClass":"pl-c1"},{"start":10,"end":11,"cssClass":"pl-k"},{"start":11,"end":17,"cssClass":"pl-en"},{"start":18,"end":22,"cssClass":"pl-smi"},{"start":23,"end":32,"cssClass":"pl-smi"}],[{"start":4,"end":8,"cssClass":"pl-en"},{"start":10,"end":16,"cssClass":"pl-s"},{"start":10,"end":11,"cssClass":"pl-pds"},{"start":15,"end":16,"cssClass":"pl-pds"},{"start":17,"end":24,"cssClass":"pl-en"},{"start":25,"end":27,"cssClass":"pl-smi"}],[],[{"start":11,"end":12,"cssClass":"pl-k"},{"start":13,"end":17,"cssClass":"pl-smi"},{"start":18,"end":27,"cssClass":"pl-en"},{"start":28,"end":30,"cssClass":"pl-smi"}],[],[{"start":10,"end":11,"cssClass":"pl-k"},{"start":12,"end":18,"cssClass":"pl-smi"}],[{"start":12,"end":13,"cssClass":"pl-k"},{"start":14,"end":19,"cssClass":"pl-en"},{"start":20,"end":24,"cssClass":"pl-smi"},{"start":25,"end":29,"cssClass":"pl-smi"},{"start":30,"end":34,"cssClass":"pl-smi"},{"start":35,"end":38,"cssClass":"pl-smi"}],[{"start":17,"end":18,"cssClass":"pl-k"},{"start":19,"end":24,"cssClass":"pl-en"},{"start":25,"end":33,"cssClass":"pl-smi"},{"start":34,"end":38,"cssClass":"pl-smi"},{"start":39,"end":43,"cssClass":"pl-smi"},{"start":44,"end":47,"cssClass":"pl-smi"}],[],[{"start":6,"end":7,"cssClass":"pl-k"},{"start":8,"end":9,"cssClass":"pl-c1"}],[{"start":4,"end":53,"cssClass":"pl-c"},{"start":4,"end":5,"cssClass":"pl-c"}],[{"start":4,"end":6,"cssClass":"pl-k"},{"start":7,"end":21,"cssClass":"pl-smi"},{"start":22,"end":24,"cssClass":"pl-k"},{"start":25,"end":26,"cssClass":"pl-c1"}],[{"start":12,"end":13,"cssClass":"pl-k"},{"start":14,"end":18,"cssClass":"pl-smi"},{"start":19,"end":20,"cssClass":"pl-smi"},{"start":21,"end":22,"cssClass":"pl-k"},{"start":23,"end":24,"cssClass":"pl-c1"},{"start":24,"end":25,"cssClass":"pl-k"},{"start":25,"end":26,"cssClass":"pl-smi"}],[{"start":13,"end":14,"cssClass":"pl-k"},{"start":15,"end":18,"cssClass":"pl-en"},{"start":20,"end":21,"cssClass":"pl-k"},{"start":21,"end":23,"cssClass":"pl-smi"},{"start":23,"end":25,"cssClass":"pl-k"},{"start":25,"end":26,"cssClass":"pl-c1"},{"start":27,"end":28,"cssClass":"pl-k"},{"start":29,"end":32,"cssClass":"pl-smi"},{"start":37,"end":45,"cssClass":"pl-c"},{"start":37,"end":39,"cssClass":"pl-c"},{"start":40,"end":45,"cssClass":"pl-en"}],[{"start":22,"end":23,"cssClass":"pl-k"},{"start":24,"end":25,"cssClass":"pl-c1"},{"start":25,"end":26,"cssClass":"pl-k"},{"start":26,"end":30,"cssClass":"pl-en"},{"start":31,"end":35,"cssClass":"pl-smi"},{"start":36,"end":37,"cssClass":"pl-c1"}],[{"start":8,"end":21,"cssClass":"pl-en"},{"start":22,"end":33,"cssClass":"pl-smi"},{"start":35,"end":36,"cssClass":"pl-k"}],[{"start":18,"end":19,"cssClass":"pl-k"},{"start":20,"end":55,"cssClass":"pl-en"},{"start":57,"end":61,"cssClass":"pl-smi"},{"start":64,"end":65,"cssClass":"pl-smi"},{"start":66,"end":67,"cssClass":"pl-k"},{"start":68,"end":75,"cssClass":"pl-smi"},{"start":78,"end":79,"cssClass":"pl-smi"},{"start":82,"end":88,"cssClass":"pl-smi"},{"start":91,"end":94,"cssClass":"pl-smi"},{"start":97,"end":102,"cssClass":"pl-smi"},{"start":105,"end":118,"cssClass":"pl-smi"},{"start":121,"end":131,"cssClass":"pl-smi"},{"start":134,"end":141,"cssClass":"pl-smi"}],[{"start":21,"end":22,"cssClass":"pl-k"},{"start":23,"end":32,"cssClass":"pl-smi"},{"start":33,"end":41,"cssClass":"pl-smi"}],[{"start":4,"end":8,"cssClass":"pl-k"}],[{"start":8,"end":76,"cssClass":"pl-c"},{"start":8,"end":9,"cssClass":"pl-c"}],[{"start":8,"end":13,"cssClass":"pl-k"},{"start":18,"end":19,"cssClass":"pl-smi"},{"start":23,"end":29,"cssClass":"pl-smi"},{"start":42,"end":46,"cssClass":"pl-smi"},{"start":47,"end":48,"cssClass":"pl-smi"},{"start":51,"end":52,"cssClass":"pl-smi"},{"start":57,"end":70,"cssClass":"pl-smi"},{"start":71,"end":75,"cssClass":"pl-smi"},{"start":76,"end":83,"cssClass":"pl-smi"}],[{"start":12,"end":19,"cssClass":"pl-en"},{"start":21,"end":28,"cssClass":"pl-en"},{"start":29,"end":30,"cssClass":"pl-smi"},{"start":32,"end":35,"cssClass":"pl-s"},{"start":32,"end":33,"cssClass":"pl-pds"},{"start":34,"end":35,"cssClass":"pl-pds"}],[{"start":16,"end":17,"cssClass":"pl-k"},{"start":18,"end":22,"cssClass":"pl-smi"},{"start":23,"end":24,"cssClass":"pl-smi"},{"start":25,"end":26,"cssClass":"pl-k"},{"start":27,"end":28,"cssClass":"pl-c1"},{"start":28,"end":29,"cssClass":"pl-k"},{"start":29,"end":30,"cssClass":"pl-smi"}],[{"start":17,"end":18,"cssClass":"pl-k"},{"start":19,"end":22,"cssClass":"pl-en"},{"start":24,"end":25,"cssClass":"pl-k"},{"start":25,"end":27,"cssClass":"pl-smi"},{"start":27,"end":29,"cssClass":"pl-k"},{"start":29,"end":30,"cssClass":"pl-c1"},{"start":31,"end":32,"cssClass":"pl-k"},{"start":33,"end":36,"cssClass":"pl-smi"},{"start":41,"end":49,"cssClass":"pl-c"},{"start":41,"end":43,"cssClass":"pl-c"},{"start":44,"end":49,"cssClass":"pl-en"}],[{"start":26,"end":27,"cssClass":"pl-k"},{"start":28,"end":29,"cssClass":"pl-c1"},{"start":29,"end":30,"cssClass":"pl-k"},{"start":30,"end":34,"cssClass":"pl-en"},{"start":35,"end":39,"cssClass":"pl-smi"},{"start":40,"end":41,"cssClass":"pl-c1"}],[{"start":12,"end":25,"cssClass":"pl-en"},{"start":26,"end":37,"cssClass":"pl-smi"},{"start":39,"end":40,"cssClass":"pl-k"}],[{"start":22,"end":23,"cssClass":"pl-k"},{"start":24,"end":59,"cssClass":"pl-en"},{"start":61,"end":65,"cssClass":"pl-smi"},{"start":68,"end":69,"cssClass":"pl-smi"},{"start":70,"end":71,"cssClass":"pl-k"},{"start":72,"end":79,"cssClass":"pl-smi"},{"start":82,"end":83,"cssClass":"pl-smi"},{"start":86,"end":92,"cssClass":"pl-smi"},{"start":95,"end":98,"cssClass":"pl-smi"},{"start":101,"end":106,"cssClass":"pl-smi"},{"start":109,"end":122,"cssClass":"pl-smi"},{"start":125,"end":135,"cssClass":"pl-smi"},{"start":138,"end":145,"cssClass":"pl-smi"}],[{"start":16,"end":17,"cssClass":"pl-k"},{"start":18,"end":25,"cssClass":"pl-smi"},{"start":26,"end":27,"cssClass":"pl-k"},{"start":28,"end":37,"cssClass":"pl-smi"},{"start":38,"end":41,"cssClass":"pl-smi"}],[{"start":21,"end":22,"cssClass":"pl-k"},{"start":23,"end":35,"cssClass":"pl-smi"},{"start":36,"end":37,"cssClass":"pl-k"},{"start":38,"end":47,"cssClass":"pl-smi"},{"start":48,"end":56,"cssClass":"pl-smi"}],[{"start":14,"end":15,"cssClass":"pl-k"},{"start":16,"end":17,"cssClass":"pl-smi"},{"start":17,"end":18,"cssClass":"pl-k"},{"start":18,"end":19,"cssClass":"pl-c1"}],[{"start":20,"end":21,"cssClass":"pl-k"},{"start":22,"end":25,"cssClass":"pl-smi"}],[{"start":25,"end":26,"cssClass":"pl-k"},{"start":27,"end":35,"cssClass":"pl-smi"}],[{"start":8,"end":11,"cssClass":"pl-k"}],[{"start":20,"end":21,"cssClass":"pl-k"},{"start":22,"end":26,"cssClass":"pl-smi"},{"start":27,"end":34,"cssClass":"pl-smi"}],[{"start":19,"end":21,"cssClass":"pl-smi"},{"start":23,"end":24,"cssClass":"pl-k"},{"start":26,"end":30,"cssClass":"pl-en"},{"start":31,"end":35,"cssClass":"pl-smi"},{"start":36,"end":37,"cssClass":"pl-smi"},{"start":37,"end":38,"cssClass":"pl-k"},{"start":38,"end":39,"cssClass":"pl-c1"},{"start":41,"end":42,"cssClass":"pl-smi"},{"start":42,"end":43,"cssClass":"pl-k"},{"start":43,"end":44,"cssClass":"pl-c1"},{"start":47,"end":48,"cssClass":"pl-smi"},{"start":48,"end":49,"cssClass":"pl-k"},{"start":49,"end":50,"cssClass":"pl-c1"}],[{"start":4,"end":7,"cssClass":"pl-k"}],[],[{"start":12,"end":14,"cssClass":"pl-smi"},{"start":16,"end":17,"cssClass":"pl-k"},{"start":18,"end":30,"cssClass":"pl-smi"}],[{"start":4,"end":11,"cssClass":"pl-en"},{"start":12,"end":16,"cssClass":"pl-s"},{"start":12,"end":13,"cssClass":"pl-pds"},{"start":13,"end":15,"cssClass":"pl-cce"},{"start":15,"end":16,"cssClass":"pl-pds"}],[{"start":0,"end":3,"cssClass":"pl-k"}],[],[{"start":9,"end":10,"cssClass":"pl-k"},{"start":11,"end":15,"cssClass":"pl-en"},{"start":16,"end":20,"cssClass":"pl-smi"},{"start":21,"end":22,"cssClass":"pl-smi"}]],"csv":null,"csvError":null,"dependabotInfo":{"showConfigurationBanner":null,"configFilePath":null,"networkDependabotPath":"/nicolasduchateau/multiscale-kernel-regression/network/updates","dismissConfigurationNoticePath":"/settings/dismiss-notice/dependabot_configuration_notice","configurationNoticeDismissed":false,"repoAlertsPath":"/nicolasduchateau/multiscale-kernel-regression/security/dependabot","repoSecurityAndAnalysisPath":"/nicolasduchateau/multiscale-kernel-regression/settings/security_analysis","repoOwnerIsOrg":false,"currentUserCanAdminRepo":true},"displayName":"interpolateMulti.m","displayUrl":"https://github.com/nicolasduchateau/multiscale-kernel-regression/blob/main/interpolateMulti.m?raw=true","headerInfo":{"blobSize":"4.24 KB","deleteInfo":{"deleteTooltip":"Delete this file"},"editInfo":{"editTooltip":"Edit this file"},"ghDesktopPath":"https://desktop.github.com","gitLfsPath":null,"onBranch":true,"shortPath":"e5bc914","siteNavLoginPath":"/login?return_to=https%3A%2F%2Fgithub.com%2Fnicolasduchateau%2Fmultiscale-kernel-regression%2Fblob%2Fmain%2FinterpolateMulti.m","isCSV":false,"isRichtext":false,"toc":null,"lineInfo":{"truncatedLoc":"129","truncatedSloc":"113"},"mode":"executable file"},"image":false,"isCodeownersFile":null,"isPlain":false,"isValidLegacyIssueTemplate":false,"issueTemplateHelpUrl":"https://docs.github.com/articles/about-issue-and-pull-request-templates","issueTemplate":null,"discussionTemplate":null,"language":"MATLAB","languageID":225,"large":false,"loggedIn":true,"newDiscussionPath":"/nicolasduchateau/multiscale-kernel-regression/discussions/new","newIssuePath":"/nicolasduchateau/multiscale-kernel-regression/issues/new","planSupportInfo":{"repoIsFork":null,"repoOwnedByCurrentUser":null,"requestFullPath":"/nicolasduchateau/multiscale-kernel-regression/blob/main/interpolateMulti.m","showFreeOrgGatedFeatureMessage":null,"showPlanSupportBanner":null,"upgradeDataAttributes":null,"upgradePath":null},"publishBannersInfo":{"dismissActionNoticePath":"/settings/dismiss-notice/publish_action_from_dockerfile","dismissStackNoticePath":"/settings/dismiss-notice/publish_stack_from_file","releasePath":"/nicolasduchateau/multiscale-kernel-regression/releases/new?marketplace=true","showPublishActionBanner":false,"showPublishStackBanner":false},"renderImageOrRaw":false,"richText":null,"renderedFileInfo":null,"shortPath":null,"tabSize":8,"topBannersInfo":{"overridingGlobalFundingFile":false,"globalPreferredFundingPath":null,"repoOwner":"nicolasduchateau","repoName":"multiscale-kernel-regression","showInvalidCitationWarning":false,"citationHelpUrl":"https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/about-citation-files","showDependabotConfigurationBanner":null,"actionsOnboardingTip":null},"truncated":false,"viewable":true,"workflowRedirectUrl":null,"symbols":{"timedOut":false,"notAnalyzed":true,"symbols":[]}},"copilotInfo":{"documentationUrl":"https://docs.github.com/copilot/overview-of-github-copilot/about-github-copilot-for-business","notices":{"codeViewPopover":{"dismissed":false,"dismissPath":"/settings/dismiss-notice/code_view_copilot_popover"}},"userAccess":{"accessAllowed":false,"hasSubscriptionEnded":false,"orgHasCFBAccess":false,"userHasCFIAccess":false,"userHasOrgs":true,"userIsOrgAdmin":false,"userIsOrgMember":false,"business":null,"featureRequestInfo":null}},"csrf_tokens":{"/nicolasduchateau/multiscale-kernel-regression/branches":{"post":"5_rHSZICmEeoXJcAjpdRFi_JNR7wqn9bAK2sg1m9vcO3jDXKT7IQayio5OjLkni-Z2Ctw3bqqzcMmRtzAeMwIA"},"/repos/preferences":{"post":"BUbGEVR0EEePcPDr7rHxzOQljupYWPVIOFDSyNllVS6kl8M7jJohDXQ_itCdD-tEVvUYf9BHzrFrb2HZ59AqTw"}}},"title":"multiscale-kernel-regression/interpolateMulti.m at main · nicolasduchateau/multiscale-kernel-regression"}