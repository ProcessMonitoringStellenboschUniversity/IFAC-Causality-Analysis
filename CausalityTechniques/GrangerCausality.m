classdef GrangerCausality<handle
    %GRANGERCAUSALITY class for calculation of Granger Causality
    
    properties
        IDs
        ConnMatrix
        NewCM
        PVals
        Parameters
        
    end
    
    methods
        
        function this = GrangerCausality()
            %Granger Causality Class constructor
            
        end
        
        function this = calculate(this,data,names)
            
            [cm, pval,~,k,~] = this.grangerCaus(data);
            this.refine(cm,pval,0.05)
            % Package results
            this.ConnMatrix = cm;
            this.PVals = pval;
            %             this.Parameters.KMax = kMax;
            this.Parameters.K = k;
            this.IDs = names;
            
        end
        
        function [cm, pval,fStat,K,aic] = grangerCaus(this,X,varargin)
            %GRANGERCAUS Granger causality
            
            %OUTPUTS: cm: MxM causality matrix. Entry i,j indicates the
            %causal influence from variable i to variable j
            
            [N,M] = size(X);
            
            if nargin>2
                % if model order is specified
                K = varargin{1};
                aic = [];
            else
                % Model order estimation
                kMax = 40;
                aic = zeros(kMax,1);
                SIGMAfull = cell(kMax,1);
                E = cell(kMax,1);
                for k = 1:kMax
                    % Get residuals covariance matrix
                    [~,SIGMAfull{k},E{k}] = this.dataToVar(X,k);
                    % Using the AIC from Duan (2014) (Eq 4.18)
                    fprintf('Model order = %d\n',k);
                    % working out aic with just sigma from X to Y
                    aic(k) = (N-k)*log(det(SIGMAfull{k}(1))) + (N-k)*(2*k*M^2)/((N-k)-(k*M^2)-1);%Toolbox uses a corrected version
                    
                end
                % Select k with lowest AIC value:
                [~,K] = min(aic);
                fprintf('Chosen model order: K = %d\n',K);
                %                 figure
                %                 plot(1:kMax,aic,'bo-')
                %                 title('Model Order Esitmation: Akaike Information Criteria')
                %                 xlabel('Model Order')
                %                 ylabel('AIC')
                
            end
            % Full regression
            [~,SIGMAfull,E] = this.dataToVar(X,K);
            rssF = sum(E.^2);
            
            % Reduced regression and F claculation:
            vars = 1:M;
            cm = zeros(M,M);
            fStat = zeros(M,M);
            pval = zeros(M,M);
            lSigFull = log(diag(SIGMAfull));                       
            for r = 1:M
                % Regression including all variables except variable No. r
                omitVar = vars(vars~=r); % omit variable No. r
                [~, sigma,eR] = this.dataToVar(X(:,omitVar),K);
                lSig = log(diag(sigma));
                rssR = sum(eR.^2);
                for ind = 1:M-1
                    c = omitVar(ind);
                    %influence omitted variable has on variable no. c
                    cm(r,c) = lSig(ind)-lSigFull(c);
                    fStat(r,c) = (rssR(ind) - rssF(c))*(N-2*K-1)/(rssF(c)*K);
                    pval(r,c) = 1-fcdf(fStat(r,c),K,(N-2*K-1));
                end
            end
        end
        
        function [cm, pVal,fStat,K] = grangerCausWrap(this,X,varargin)
            %GRANGERCAUS Pairwise granger causality (nonconditional)
            
            [~,M] = size(X);
            cm = zeros(M,M);
            pVal = zeros(M,M);
            fStat = zeros(M,M);
            for r = 1:M
                for c = 1:M
                    if c == r
                        % Skipping diagonals
                    else
                        if nargin>2
                            % if model order is specified
                            K = varargin{1};
                            [gccm, pval,fstat,~,~] = this.grangerCaus([X(:,r) X(:,c)],K);
                        else
                            [gccm, pval,fstat,K(r,c),~] = this.grangerCaus([X(:,r) X(:,c)]);
                        end
                        cm(r,c) = gccm(1,2);
                        %                         cm(c,r) = gccm(1,2);
                        pVal(r,c) = pval(1,2);
                        %                         pVal(c,r) = pval(1,2);
                        fStat(r,c) = fstat(1,2);
                        %                         fStat(c,r) = fstat(1,2);
                    end
                end
            end
        end
                
    end
    
    methods (Static)
        
        function [ aic, aicOrder ] = orderEstimation(X,kMax)
            %ORDERESTIMATION Estimaiton of model order using Akaike information
            %criteria
            
            [N,M] = size(X);
            aic = zeros(kMax,1);
            for k = 1:kMax
                % Get residuals covariance matrix
                [~,SIGMAfull,~] = dataToVar(X,k);
                % Using the AIC from Duan (2014) (Eq 4.18)
                fprintf('Model order = %d\n',k);
                aic(k) = log(det(SIGMAfull)) + (2*k*M^2)/(N-k);%Toolbox uses a corrected version
            end
            [~,aicOrder] = min(aic);
            figure
            plot(1:kMax,aic,'bo-')
            title('Model Order Esitmation: Akaike Information Criteria')
            xlabel('Model Order')
            ylabel('AIC')
        end
        
        function [A,SIGMA,E] = dataToVar(X,K)
            %
            %   Detailed explanation goes here
            
            [N,M] = size(X);
            K1 = K+1;
            X0 = X(K1:end,:);
            XL = zeros(N-K,M,K);
            for k = 1:K
                XL(:,:,k) = X(K1-k:N-k,:);
            end
            XL = reshape(XL,N-K,M*K);
            A = X0'/XL';
            E = X0'-A*XL';
            SIGMA = (E*E')/(N-K-1);
            E = E';
            A = reshape(A,M,M,K);
        end
    end
end

