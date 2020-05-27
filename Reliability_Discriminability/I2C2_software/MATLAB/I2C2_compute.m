function [lambda, trKu, trKx, W] = I2C2_compute(y, varargin)
%computeI2C2 computes the I2C2 value for the given matrix y of data
%   Detailed explanation goes here
% y: An n by p data matrix containing functional responses. Each row contains measurements 
%    from a function for one observation at a set of grid points, and each column contains measurements
%    of all functions at a particular grid point. 
%    The rows are organized by subjects and then visits, EX) (Y11, Y12, Y21, Y22, ... , YI1 , YI2)
% I: Number of subjects  
% J: Number of repetitions
% id: Vector (column) of IDs, EX)c(1, 1, 2, 2, 3, 3, 4, 4, ... , I, I)
% visit: Vector of visits (column), EX) (1, 2, 1, 2, 1, 2, ... , 1, 2)
% T: Vector of distance (in time) of each visit from baseline visit.
%	   If the length of time between visits is different for the subjects in the dataset
%    then match visits according to their distance in time as given by T. 
%    If T == NULL, match observations from different clusters by visit number
% twoway: a logical argument indicating whether a oneway or a twoway functional ANOVA (analysis
%    of variance) decomposition is more appropriate for the problem. "twoway = TRUE" will carry
%    out twoway ANOVA and remove the visit specific mean
% demean: if TRUE, include the demean step and output the demeaned dataset
% symmetric: if FALSE then the function uses the method of moments estimator formula; 
%    if TRUE, pairwise symmetric sum formula, 
%    default is FALSE
% trun: if TRUE, set negative I2C2 to zero
% 
%Output
% 
%     The output of the function is a list that contains the following elements. 
%
% lambda:       estimated I2C2
% Kx:           the trace of between-cluster variance operator
% Ku:           the trace of within-cluster variance operator
% W:     if demean == TRUE, output the demeaned dataset, otherwise
% the original dataset
%
% Author(s): 
%
%    Adina Crainiceanu, Haochang Shou, Ani Eloyan, Seonjoo Lee, Vadim Zipunnikov, Mary Beth Nebel, 
%    Brian Caffo, Martin Lindquist, Ciprian M. Crainiceanu     
%
%References:
%
%    Haochang Shou, Ani Eloyan, Seonjoo Lee, Vadim Zipunnikov, Mary Beth Nebel, 
%    Brian Caffo, Martin Lindquist, Ciprian M. Crainiceanu  (2012) The image intra 
%    class correlation coefficient for replication studies.
%

%check for optional params
if nargin < 5
    error('at least 5 arguments, (y, "I", I, "J", J), or (y, "id", id, "visit", visit) required');
end

% Default settings
demean  = true;
twoway = true;
symmetric = false;
trun = false;
I=0;
J=0;
id = [];
visit = [];
T = [];

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

if (demean)
    %compute mean per column and subtract it from each element
    mu = mean(y); %same as mean(y,1), produces a row with the means per column
    mu_matrix = repmat(mu, n, 1); % repeat the row n times, so we get a matrix
    resd = y - mu_matrix; % subtract the means from the original matrix
    % y = y-mu;
    
    % If twoway functional ANOVA is needed ("twoway==TRUE"),  the visit specific mean function
    % and the deviation from the overall mean to visit specific mean functions are also
    % computed. 
    if( twoway) 
	       if isempty(T)  
	     	   T = visit;
           end
		   eta = zeros(length(unique(T)), p);
           for j = unique(T)' %make sure we iterate through a row, not column, otherwise j will take the entire vector value
			   if (sum( T == j ) ) == 0 
                   continue;
               end
			   if sum( T ==j ) == 1   
                   size(T)
                   size(eta)
                   size(y)
                   size(y(T==j,:))
                   size(mu)
                   size(eta( unique(T) == j,:))
			  	   %TODO debug here the sizes do not agree
                   T
                   unique(T)
                   j
                   eta( unique(T) == j,:) = y(T == j, : ) - mu;
               else
                   eta( unique(T) == j , :) = mean( y( T == j, :), 1) - mu;
               end
           end

           % Calculate residuals by subtracting visit-specific mean from original functions for 
        % 'twoway == TRUE', or subtracting overall mean function for 'twoway == FALSE'. 

		   for j = unique(T)'   
			   if sum( T == j ) == 0 
                   continue;
               end
               mu_eta = mu + eta( unique(T) == j ,:);
               mu_eta_matrix = repmat(mu_eta,size(y(T==j,:),1),1); %same as mu_eta(ones(y(T==j,:),1),1),:)
		 	   resd(T == j,:) = y( T == j, :) - mu_eta_matrix;
           end				
      end
    
    
    W = resd;
else
    W = y;
end

%reset the ids to be from 1 to I 
[~, id] = ismember(id, unique(id, 'stable'));

%compute k2 - the sum for all i of J_i^2 (sum of squares of number of
%visits for each patient)
n_I0_all = tabulate(id); 
%tabulate gives value, number of times for the value in the vector, 
%percentage in the vector

n_I0 = n_I0_all(:,2); %take only the number of visits for each patient, 
%size(n_I0)
% size n_I0 is (I,1)

k2=sum(n_I0.^2);

%compute W dot dot - population average
Wdd = mean(W,1); %same as mean(W,1), produces a row with the means per column

%compute subject-specific means Si (I x p) matrix
Si = grpstats(W, id, 'sum'); %compute the sum by subject 
Si_mean = grpstats(W, id, 'mean'); %compute the mean by subject 
%W(:,1)
%Si(:, 1)
%id
% return;
%if symmetric is FALSE, use the method of moments estimator formula from
%the manuscript, with a correction for bias by Seoungee

if (~symmetric)
 
    %duplicate the rows so we still get a nxp matrix from the per-subject means
    Wi_matrix = repelem(Si_mean,n_I0,1); 
    %repelem (matrix, vector, 1) repeats each row 
    % the number of times specified in the vector position for that row
    
    %duplicate the rows so we still get a nxp matrix from the population means
    Wdd_matrix = kron(Wdd,ones(n,1));

    %compute trace Kw and Ku
    trKw = sum(sum((W-Wdd_matrix).^2))/(n-1);
    trKu = sum(sum((W-Wi_matrix).^2))/(n-I);

    %compute Kx, but the the bias correction factor needed when symmetric is
    %false
     trKx = (trKw-trKu)/(1+(1-k2/n)/(n-1));
    %trKx = (trKw-trKu);
 
else
    n_I0_id = repelem(n_I0,n_I0,1);
    %repelem (matrix, vector, 1) repeats each row 
    % the number of times specified in the vector position for that row
    %TODO - this repmat to make a matrix out of the column is waste of
    %memory
    % a loop will probably be much faster to multiply each column with the
    % n_10_id
    n_I0_id_matrix = repmat(n_I0_id,1,p);
    
    %compute trace Kw and Ku
    trKu = (sum(sum(W.^2 .*n_I0_id_matrix))-sum(sum(Si.^2)))/(k2-n);
    %size(trKu)
     trKw = (sum(sum(W.^2))*n - sum(sum((n*Wdd).^2)) - trKu *(k2-n))/(n^2-k2);
     %size(trKw)
     %compute trace Kx
     trKx = trKw - trKu;
    %size(trKx)
 end
  

lambda = trKx/(trKx+trKu); %estimated I2C2 value
if trun && lambda < 0
    lambda = 0;
end

I2C2= lambda;

;
end

