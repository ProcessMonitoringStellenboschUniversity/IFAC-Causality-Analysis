classdef TransferEntropy<handle
    %TRANSFERENTROPY class for calculation of transfer entropy
    
    properties
        
        ConnMatrix
        NewCM
        IDs
        SigThresholds
        Parameters
        
    end
    
    methods
        
        function this = TransferEntropy()
            %TRANSFERENTROPY Class constructor
            
        end
        
        function this = calculate(this,data,names)
                        
           cm = this.transEntropy(data);
           S = this.surrogateThreshMatrix(data,100,4,1,2);
           [cmRef] = this.refineCM(cm,S);
            % Package results
            this.ConnMatrix = cm;
            this.NewCM = cmRef;
            this.SigThresholds = S;
            this.IDs = names;
%             this.Parameters.KMax = kMax;
%             this.Parameters.K = k;
        end
        
        function [Txy, varargout] = transEntropy(this,data,varargin)
            %TRANSENTROPY Calculate transfer entropy. Straightforward (not
            %direct transfer entropy)
            
            [~, M] = size(data);
            nBins = 10; % Number of bins for probability density estimation
            % Parameters for transfer entropy calculation
            if nargin>2
                H = varargin{1};
                tau = H;
                K = varargin{2};
                L = varargin{3};
            else
                % Default values (from Bauer et al. 2007)
                H = 4; %prediction horizon.
                tau = H;% If possible h and tau should be equal to dead time
                K = 1;% Embedding dimension for x
                L = 2;% Embedding dimension for y
            end
            
            % Scaling
            [data, ~, ~] = zscore(data);
            
            % Loop for every variable
            Hyy = zeros(M,1);
            Hyyx = zeros(M,M);
            Txy = zeros(M,M);
            for c = 1:M
                % Get conditional shannon entropy for Y with past Y
                Hyy(c) = this.getShannonY(data(:,c),K,H,tau,nBins);
                % Inner loop for every pair of variables
                for r = 1:M
                    if r==c
                        % Skipping diagonals
                    else
                        % Get conditional shannon entropy for Y with past Y and X
                        Hyyx(r,c) = this.getShannonX(data(:,c),data(:,r),K,L,H,tau,nBins);
                        Txy(r,c) = Hyy(c) - Hyyx(r,c);% Txy(1,2) should indicate x1-->x2                        
                    end
                end
            end

            if nargout>1
                varargout{1} = Hyy;
            end
            if nargout>2
                varargout{2} = Hyyx;
            end
            
        end
        
        function [Txy, varargout] = transEntropySinglePair(this,X,Y,varargin)
            %TRANSENTROPY Calculate transfer entropy. Straightforward (not
            %direct transfer entropy)
            nBins = 10; % Number of bins for probability density estimation
            % Parameters for transfer entropy calculation
            if nargin>3
                H = varargin{1};
                tau = varargin{2};
                K = varargin{3};
                L = varargin{4};
            else
                % Default values (from Bauer et al. 2007)
                H = 4; %prediction horizon.
                tau = H;% If possible h and tau should be equal to dead time
                K = 1;% Embedding dimension for x
                L = 2;% Embedding dimension for y
            end
            
            % Scaling
            [X, ~, ~] = zscore(X);
            [Y, ~, ~] = zscore(Y);
            
            % Get conditional shannon entropy for Y with past Y
            Hyy = this.getShannonY(Y,K,H,tau,nBins);
            % Inner loop for every pair of variables
            
            % Get conditional shannon entropy for Y with past Y and X
            Hyyx = this.getShannonX(Y,X,K,L,H,tau,nBins);
            Txy = Hyy - Hyyx;

            if nargout>1
                varargout{1} = Hyy;
            end
            if nargout>2
                varargout{2} = Hyyx;
            end
            
        end
        
        function [Hyy] = getShannonY(this,Y,K,H,tau,nBins)
            %GETSHANNONY Conditional Shannon entropy for Y given past Y
            %Used for Transfer Entropy calculation.
            %NOTE: using notation where X is input and Y is output, I.e.
            %to calculate TE from X-->Y
            
            N = length(Y);
            Y = zscore(Y);
            yH = Y(1:end-H);
            % Creating Embedded matrix
            yEmb = zeros(N-H,K);
            for k = 1:K
                yEmb(H+k*tau:end,k) = yH(1:end-k*tau-H+1);
            end
            %Calculate probabilities
            yH = yH(K*tau+H:end);
            yi = yEmb(K*tau+H:end,1:K);
            yyi = [yH yi];
            if K<6
                [Pyi,~] = this.kpdf(yi,nBins);
                [Pyyi, gridY] = this.kpdf(yyi,nBins);
                % Calculate Shannon entropy for Y, Hyy
                hyy = Pyyi.*log2(Pyyi./repmat(Pyi,nBins,1));
                hyy = sum(hyy);
            else
                [Pyi,~] = this.kpdf(yi,nBins,'first');
                [Pyyi, gridY] = this.kpdf(yyi,nBins,'second');
                
                %                 div = 0:length(Pyi):length(Pyyi);
                
                hyy = 0;
                int = 100;
                div = 0:int:length(Pyyi.p)-int;
                div2 = mod(div,length(Pyi.p));
                %                 count = 0;
                parfor d = 1:length(div)
                    %                     hyy = hyy + sum(Pyyi(div(lp-1)+1:div(lp)).*log2(Pyyi(div(lp-1)+1:div(lp))./Pyi));
                    %                     if d == 1
                    %                         count = 0;
                    %                     end
                    %                     if count == length(Pyi.p)
                    %                         count = 0;
                    %                     end
                    hyy = hyy + sum(Pyyi.p(div(d)+1:div(d)+int,1).*log2(...
                        Pyyi.p(div(d)+1:div(d)+int,1)./Pyi.p(div2(d)+1:div2(d)+int,1)));
                    %                    count = count+int;
                end
            end
            
            
            % Sum over all amplitude bins
            %             hyy = sum(hyy);
            Hyy = -hyy*(gridY^(K+1));
        end
        
        function [Hyyx] = getShannonX(this,Y,X,K,L,H,tau,nBins)
            %GETSHANNONY Conditional Shannon entropy for Y given past Y & X
            %Used for Transfer Entropy calculation.
            %NOTE: using notation where X is input and Y is output, I.e.
            %to calculate TE from X-->Y
            N = length(Y);
            
            yH = Y(1:end-H);
            % Creating Embedded matrices
            yEmb = zeros(N-H,K);
            for k = 1:K
                yEmb(H+k*tau:end,k) = yH(1:end-k*tau-H+1);
            end
            xH = X(1:end-H);
            xEmb = zeros(N-H,L);
            for l = 1:L
                %Only need one embedded matrix, unless tau is different for
                %x and y
                xEmb(H+l*tau:end,l) = xH(1:end-l*tau-H+1);
            end
            
            %Loop for every pair of variables
            %             Hyyx = zeros(M,M);
            %             for r = 1:M
            %                 for c = 1:M
            %                     xi = squeeze(xEmb(L*tau+H:end,c,1:L));
            %                     yx = [squeeze(xEmb(L*tau+H:end,r,1:K)) xi];
            xi = xEmb(L*tau+H:end,1:L);
            yx = [yEmb(L*tau+H:end,1:K) xi];
            yyx =[yH(L*tau+H:end) yx];
            [Pyx,~] = this.kpdf(yx,nBins);
            [Pyyx, gridYX] = this.kpdf(yyx,nBins);
            % Calculate Shannon entropy of Y given X and Y
            if (K+L) ==1
                hyyx = Pyyx.*log2(Pyyx./repmat(Pyx,nBins,1));
            else
                hyyx = Pyyx.*log2(Pyyx./repmat(Pyx,nBins,1));
            end
            % Sum over all amplitude bins
            hyyx = sum(hyyx);
            Hyyx= -hyyx*(gridYX^(L+1));
            %                 end
            %             end
        end
        
       
        
    end
    
    methods (Static)
        
        function [p, gridSize] = kpdf(X, nBins,varargin)
            %KPDF Probability density function using kernel estimation
            %method
            if isempty(X)
                p = 1;
                gridSize = 1;
            else
                [N,M] = size(X);
                % rule of thumb bandwidth suggested by Bowman and Azzalini (1997) p.31
                % h=median(abs(x-repmat(median(x),N,1)))/0.6745*(1/N)^(1/6);
                % rule of thumb from Duan (2014)
                c = (4/3)^(1/5);
                h = c*std(X)*N^(-1/(M+4));
                xmax = max(X)+3*h;
                for m = 1:M;
                    % Doing this to ignore zeros in matrices
                    xmin(m) = min(X(X(:,m)~=0,m))-3*h(m);
                end
                
                xpts = cell(1,M);
                for m = 1:M
                    xpts{m}=linspace(xmin(m),xmax(m),nBins);
                end
                gridSize = xpts{1}(2) - xpts{1}(1);
                
                xGrid = cell(1,M);
                [xGrid{:}] = ndgrid(xpts{:});
                %To make grid equivalent to that obtained from meshgrid, first two outputs
                %have to be permuted
                if M>1
                    permDims = 1:1:M;
                    permDims(1) = 2; permDims(2) = 1;
                    xGrid{1} = permute(xGrid{1},permDims);
                    xGrid{2} = permute(xGrid{2},permDims);
                end
                
                %Flattened grid
                Xg =zeros(nBins^M,M);
                for m = 1:M
                    Xg(:,m) = xGrid{m}(:);
                end
                
                %                 if nBins^M<100000
                % initialising p
                p = zeros(nBins^M,1);
                
                %Get current parallel pool. If no pool it creates its own
%                 pl = gcp;
                
                % calculating pdf
                parfor nn = 1:nBins^M
                    % Not too sure why this re-evaulation of h is necessary
                    hNew = h;%min([min([Xg(nn,:)-xmin;h]);min([xmax-Xg(nn,:);h])]);
                    if prod(hNew)<eps(1/M)
                        %         pdf(k)=0; %Don't need to assign zero since it is already
                        %         initialised as zero
                    else
                        v=(repmat(Xg(nn,:),N,1)-X)./repmat(hNew,N,1);
                        % Gaussian Kernel Function
                        K = exp(-0.5*sum(v.^2,2))/((2*pi)^(0.5*M));
                        p(nn) = sum(K)/(N*prod(hNew));
                        
                    end
                end
                
                % if extra argument is passed it means that we want to
                % wrtie directly to file rather than keeping large
                % probability variable in memory
                if nargin>2
                    classdir= fileparts(mfilename('fullpath'));
                    tempvar = fullfile(classdir,'temp',['p' varargin{1} '.mat']);
                    save(tempvar,'p','-v7.3')
                    p = matfile(tempvar);
                end
                %                 end
            end
        end
  
        
    end
    
end

