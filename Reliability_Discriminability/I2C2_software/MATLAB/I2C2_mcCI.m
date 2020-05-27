%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   This file contains MATLAB codes to make inference for the estimated I2C2
%%%   I2C2_mcCI provides the confidence interval by bootstrapping subjects
%%%   I2C2_mcCI calls I2C2_rsample to conduct the bootstrap
%%%   I2C2_NullDist conducts permutation by calling I2C2_rpermute to obtain
%%%   the null distribution of I2C2
%%%         MATLAB version by Adina Crainiceanu
%%%         after R version by Haochang, Senjoo and Ani
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [lambda, CI] = I2C2_mcCI(y, varargin)


% I2C2_mcCI          
%
% Compute the Confidence Interval
%
% Description:
%
%     ‘I2C2_mcCI’ is used for computing the confidence interval of I2C2 via multicore computing.
%
%Sample usage:
%
%	I2C2_mcCI(y, 'I',I, 'J', J, 'id', [], 'visit', [], 'p', [],
%		      'twoway', TRUE, 'demean', TRUE, 'T', [], 'symmetric', FALSE, 'trun', FALSE, %trace calculation parameters
%    		  'rseed', 1234, 'R', 100, 'ncores', 1, 'ci', 0.95, % Bootstrap Arguments
%              ...)
%
%Arguments:
% % y: Dataset (Y11, Y12, ... , YnJn)' => (IJ)-by-p matrix for balanced case, EX) (Y11, Y12, Y21, Y22, ... , YI1 , YI2)
% % I: Number of subjects
% % J: Number of repetitions
% % id: Vector of IDs, EX) c(1, 1, 2, 2, 3, 3, 4, 4, ... , I, I)
% % visit: Vector of visits, EX) (1, 2, 1, 2, 1, 2, ... , 1, 2)
% % p: dimension of oberved vectors Yij, EX) Number of voxels %TODO:
% currently not used as param, assumed to be same as nb of columns in y
% % twoway, demean, T symmetric, trun: trace calcuation parameters. See 'I2C2'
% % rseed: Seed number
% % R: The bootstrap repetition size
% % ncores: Number of Cores
% % ci: 100*ci% The level of the Confidence Interval
%
%Author(s):
%
%    Adina Crainiceanu,Haochang Shou, Ani Eloyan, Seonjoo Lee, Vadim Zipunnikov, Mary Beth Nebel,
%     Brian Caffo, Martin Lindquist, Ciprian M. Crainiceanu
%
%References:
%
%     Haochang Shou, Ani Eloyan, Seonjoo Lee, Vadim Zipunnikov, Adina Crainiceanu, Mary Beth Nebel,
%     Brian Caffo, Martin Lindquist, Ciprian M. Crainiceanu  (2012) The image intra
%     class correlation coe?cient for replication studies.
%
%See Also:
%
%     'I2C2_compute', 'I2C2_mcNulldist'
%
%Examples:
%
%
%%     1. Either (id,visit) or (I,J) is needed for identifying clusters.
%%         If both are missing, return an error message "not enough information";
%%         if only (I,J) is provided, the function generates (id,visit);
%%         if only ((id,visit) is provided, the function generates (I,J);
%%      Note:  currently the function only deals with balanced data, i.e, each
%%             subject is measured at the same J visits.
%%
%%     2. If the number of grid points p is missing, the function generates it.



%%%%%%%%%%%%

%check for optional params
if nargin < 5
    error('at least 5 arguments, (y, "I", I, "J", J), or (y, "id", id, "visit", visit) required');
end

% Default settings
%I=NULL, J=NULL, id = NULL, visit = NULL,  p = NULL,
%		              twoway = TRUE, demean = TRUE, T = NULL, symmetric = FALSE, trun = FALSE,
%                      rseed = 1234, R = 100, ncores = 4, ci = 0.95

I = 0;
J = 0;
id = [];
visit = [];
twoway = true;
demean  = true;
T = [];
symmetric = false;
trun = false;
rseed = 1234;
R = 100;
ncores = 1;
ci = 0.95;

% Switch to parse the varargin inputs
len = length(varargin);
% check "len" for even number
if mod(len,2) > 0
    error('Wrong arguments: must be name-value pairs.');
end
for i = 1:2:len
    switch lower(varargin{i})
        case 'i'
            I = varargin{i+1};
        case 'j'
            J = varargin{i+1};
        case 'id'
            id = varargin{i+1};
        case 'visit'
            visit = varargin{i+1};
        case 'demean'
            demean=varargin{i+1};
        case 'twoway'
            twoway=varargin{i+1};
        case 'symmetric'
            symmetric=varargin{i+1};
        case 'trun'
            trun = varargin{i+1};
        case 't'
            T = varargin{i+1};
        case 'rseed'
            rseed = varargin{i+1};
        case 'r'
            R = varargin{i+1};
        case 'ncores'
            ncores = varargin{i+1};
        case 'ci'
            ci = varargin{i+1};
        otherwise
            error(strcat('unknown parameter name ', varargin{i}));
    end
end

%either (id, visit) or (I,J) are needed. If both missing, print error
%message
if ((isempty(id) || isempty(visit)) && (I==0 || J==0))
    error('Not enough information! Please provid (I,J) or (id, visit)');
end

%if both (id, visit) and (I,J) are provided, check that they are consistent
%with each other - should be I unique values in id, each of them repeatead
%J times
if (~(isempty(id) || isempty(visit)) && ~(I==0 || J==0))
    id_tab = tabulate(id);
    if (length(id_tab) ~= I) || (sum(id_tab(:,2) ~= J) > 0)
        error('Inconsistent information from (id, visit) and (I,J)');
    end
end

%if only I,J fully provided, assume balanced design
if ((isempty(id) || isempty(visit)) && ~(I==0 || J==0))
    %create the id column 1 1 2 2 3 3 ... if J = 2
    id = repmat(1:I, J, 1);
    id = id(:);
    
    visit = repmat(1:J, I, 1);
    visit = visit';
    visit = visit(:);
end

%if only (id, visit) fully provided, compute I and J
if (~(isempty(id) || isempty(visit)) && (I==0 || J==0))
    I = length(unique(id));
    J = length(unique(visit));
end

%compute n and p the dimensions of the matrix
[n,p]=size(y);

%check that I*J = n
if I*J ~= n
    error('something wrong - I, J and nb rows in y do not match');
end

if isempty(T)
    T = visit;
end
%%	set up the number of multicores
%TODO: see how to use multiple cores in MATLAB

%	options(cores = ncores)
% 	print( paste("Number of Cores Specficed:", getOption('cores')) )

%%	Set a random seed

rng(rseed);

if ncores == 1
    lambda = arrayfun(@(x) I2C2_rsample(x),1+rseed:R+rseed);
    
end
%TODO	else {
%	lambda <- mclapply( 1:R,
%	                    myfunc <- function(s){
%	                    	                  l = I2C2_rsample(s + rseed, y = y, I = I, J = J, id = id, visit = visit, p = p,
%	                    	                                    twoway = twoway, demean = demean,
%	                    	                                    T = T, symmetric = symmetric, trun = trun
%	                    	                                    );
%	                    	                  return(l)
%	                    	                  }
%	                  )
%   }
CI = quantile( lambda, [(1 - ci)/2,ci + (1 - ci)/2] );


    function lambda = I2C2_rsample(s)
        
        
        % I2C2_rsample                 
        %
        % Compute I2C2 for a set of Randomly Sampled Subjects (with replacement)
        %
        % Description:
        %
        %     ‘I2C2_rsample’ is used for computing I2C2 of permuted data.
        %
        %Usage:
        %
        %	I2C2_rsample(s, y, I, J, id = NULL, visit = NULL, p = NULL,
        %	        	 twoway = TRUE, demean = TRUE, T = NULL, symmetric = FALSE, trun = FALSE %trace calculation parameters
        %                ...)
        %
        %Arguments:
        % % s: Random seed, default = 1
        % % y: Dataset (Y11, Y12, ... , YnJn)' => (IJ)-by-p matrix for balanced case, EX) (Y11, Y12, Y21, Y22, ... , YI1, YI2)
        % % I: Number of subjects
        % % J: Number of repetitions
        % % id: Vector of IDs, EX) c(1, 1, 2, 2, 3, 3, 4, 4, ... , I, I)
        % % visit: Vector of visits, EX) (1, 2, 1, 2, 1, 2, ... , 1, 2)
        % % p: dimension of oberved vectors Yij, EX) Number of voxels
        % % twoway, demean, T, symmetric, trun: trace calcuation parameters. See 'I2C2'
        
        %s
        rng(s);
        newI=randsample(unique(id),length(unique(id)),true);
        %pre-alocate the arrays so the size does not increase at each iteration through the next for loop
        newid=zeros(size(id));
        newX=zeros(size(y));
        newT=zeros(size(T));
        newvisit=zeros(size(visit));
        
     
        rowcounter = 1;
        for j = 1:length(unique(id))
            tmpindx=find(id==newI(j));
            nbVisits = length(tmpindx);
            newid(rowcounter:rowcounter+nbVisits-1,:) = j*ones(nbVisits,1);
            newX(rowcounter:rowcounter+nbVisits-1,:)=y(tmpindx,:);
            newT(rowcounter:rowcounter+nbVisits-1,:)=T(tmpindx);
            newvisit(rowcounter:rowcounter+nbVisits-1,:)=visit(tmpindx);
            rowcounter = rowcounter+nbVisits;
            
            %newid2=[newid2; j*ones(length(tmpindx),1)];
            %newX2=[newX2;y(tmpindx,:)];
            %newT2=[newT2;T(tmpindx)];
            %newvisit2=[newvisit2;visit(tmpindx)];
        end
        
        
        %TODO the R function has a parameter p: dimension of oberved vectors
        %Yij, EX) Number of voxels. That should be equal to nb cols in
        %newX, so I did not use the param. Check if the R fundtion works
        %when p is different (smaller) than the nb cols in newX
        
        [lambda,~,~,~]=I2C2_compute(newX, 'id',newid, 'visit',newvisit, 'J',0,'I',0, 'T',newT, 'twoway',twoway,'demean',demean,'symmetric',symmetric,'trun',trun);
        
        
    end


end