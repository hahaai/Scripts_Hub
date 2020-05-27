function lambda = I2C2_mcNulldist(y, varargin)

% I2C2_mcNulldist                   
%
% Compute Null Distribution
%
% Description:
%
%     ‘I2C2_mcNulldist’ is used for computing the null distribution of I2C2 via permutation and multicore computing.
%
%Sample usage:
%
%	I2C2_mcNulldist(y, 'I',I, 'J', J, 'id', [], 'visit', [], 'p', [],
%		      'twoway', TRUE, 'demean', TRUE, 'T', [], 'symmetric', FALSE, 'trun', FALSE, %trace calculation parameters
%    		  'rseed', 1234, 'R', 100, 'ncores', 1, % Bootstrap Arguments
%              ...)
%
%Arguments:
% % y: Dataset (Y11, Y12, ... , YnJn)' => (IJ)-by-p matrix for balanced case, EX) (Y11, Y12, Y21, Y22, ... , YI1, YI2)
% % I: Number of subjects
% % J: Number of repetitions
% % id: Vector of IDs, EX) c(1, 1, 2, 2, 3, 3, 4, 4, ... , I, I)
% % visit: Vector of visits, EX) (1, 2, 1, 2, 1, 2, ... , 1, 2)
% % p: dimension of oberved vectors Yij, EX) Number of voxels
% % twoway, demean, T, symmetric, trun: trace calcuation parameters. See 'I2C2'
% % rseed: Seed number
% % R: The bootstrap repetition size
% % ncores: Number of Cores
%
%Author(s):
%
%    Adina Crainiceanu, Haochang Shou, Ani Eloyan, Seonjoo Lee, Vadim Zipunnikov,  Mary Beth Nebel,
%     Brian Caffo, Martin Lindquist, Ciprian M. Crainiceanu
%
%References:
%
%     Haochang Shou, Ani Eloyan, Seonjoo Lee, Vadim Zipunnikov, Adina Crainiceanu, Mary Beth Nebel,
%     Brian Caffo, Martin Lindquist, Ciprian M. Crainiceanu  (2012) The image intra
%     class correlation coeffcient for replication studies.
%
%See Also:
%
%     'I2C2', 'I2C2_mcCI'
%

%check for optional params
if nargin < 5
    error('at least 5 arguments, (y, "I", I, "J", J), or (y, "id", id, "visit", visit) required');
end

% Default settings
%I = NULL, J = NULL, id = NULL, visit = NULL, p = NULL,
%			                twoway = TRUE,  demean = FALSE, T = NULL, symmetric = FALSE, trun = FALSE,
%  			                rseed = 1234, R = 500, ncores = 4

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




%TODO multiple cores
%options(cores = ncores)
%	print( paste( "Number of Cores Specficed:", getOption('cores') ) )
rng(rseed)

if ncores == 1
    
    lambda = arrayfun(@(x) I2C2_rpermute(x),1+rseed:R+rseed);
    
    %                  myfunc <- function(s){
    %                  	                   l = I2C2_rpermute(s + rseed, y = y, I = I, J = J, id = id, visit = visit, p = p,	                  	                                     twoway = twoway, demean = demean,
    %                  	                                     T = T, symmetric = symmetric, trun = trun
    %                  	                                     );
    %                  	                   return(l)
    %                  	                   }
    
    
else
    lambda = [];
    error('Multiple cores not yet supported in Matlab version of the code');
    
    %lambda <- mclapply( 1:R,
    %                    myfunc <- function(s){
    %                    	                  l = I2C2_rpermute(s + rseed, y = y, I = I, J = J, id = id, visit = visit, p = p,
    %                    	                                    twoway = twoway, demean = demean,
    %                    	                                    T = T, symmetric = symmetric, trun = trun
    %                    	                                   );
    %                    	                 return(l)
    %                    	                 }
    %                   )
end

%lambda = as.vector( unlist(lambda) )%
%return(lambda)


    function llambda = I2C2_rpermute(s)
        
        
        % I2C2_rpermute                  
        %
        % Compute I2C2 for Permuted Data
        %
        % Description:
        %
        %     ‘I2C2_rpermute’ is used for computing the I2C2 of the permuted data.
        %
        %Usage:
        %
        %	I2C2_rpermute(s, y, I, J, id = NULL, visit = NULL, p = NULL,
        %		          twoway = TRUE, demean = TRUE, T = NULL, symmetric = FALSE, trun = FALSE %trace calculation parameters
        %                 ... )
        %
        %Arguments:
        % % s: Random seed, Default = 1
        % % y: Dataset (Y11, Y12, ..., YnJn)' => (IJ)-by-p matrix for balanced case, EX) (Y11, Y12, Y21, Y22, ..., YI1, YI2)
        % % I: Number of subjects
        % % J: Number of repetitions
        % % id: Vector of IDs, EX) c(1, 1, 2, 2, 3, 3, 4, 4, ... , I, I)
        % % visit: Vector of visits, EX) (1, 2, 1, 2, 1, 2, ... , 1, 2)
        % % p: dimension of oberved vectors Yij, EX) Number of voxels
        % % twoway, demean, T, symmetric, trun: trace calcuation parameters. See 'I2C2'
        
        rng(s);
        %newI = sample( 1:nrow(y), replace = FALSE )
        newI = randsample(1:n,n, false);
        %newX = y[newI, ]
        newX = y(newI,:);
        [llambda,~,~,~]=I2C2_compute(newX, 'id',id, 'visit',visit, 'I',I, 'J',J, 'T',[], 'twoway',twoway, 'symmetric',symmetric,'demean',demean,'trun',trun);
        
    end
end
