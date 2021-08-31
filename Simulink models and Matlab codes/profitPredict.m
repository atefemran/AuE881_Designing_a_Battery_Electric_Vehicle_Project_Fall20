function [profit,demand,price] = profitPredict(attr, compAttr, betaCoeff, market)
%PROFITPREDICT predicts the profit of a product based on a linear logit model 
%
%  This linear logit model predicts [profit, demand, price] for a given vector
%  of attribute values, "attr", for a product.  By default, the profit for a
%  battery-electric vehicle is predicted. For other products, the values
%  for competition "compAttr" and linear logit coefficients "betaCoeff" must be
%  updated and provided as inputs.  The default market size is 1 million. 
%
%  The linear logit model models the demand using the following expression:
%    demand = market*exp(betaCoeff'*[price;attr(2:end)])...
%           /(exp(betaCoeff'*[price;attr(2:end)])+sum(exp(betaCoeff'*compAttr)));
%
%  This expression could have been solved directly if the price were
%  provided.  However, the first input "attr" has as its first element the
%  value for cost rather than price.  The profitPredict function therefore
%  first computes the price that maximizes profit.  This optimization
%  problem is solved analytically using the inverse Lambert W-function.
%
%  [profit,demand,price] = profitPredict(attr) returns
%     predictions for profit, demand and optimized price for a
%     battery-electric vehicle, characterized by the following
%     attr: cost in [$], range in [km], 0-100km/h-time in [s], and
%     maximum velocity in [km/h].
%     The input must be provided as a 4x1 column vector.
%
%  [profit,demand,price] = profitPredict(attr, compAttr, betaCoeff)
%     returns the outputs for a product characterized by one or more
%     competitors "compAttr" and the linear logit coefficients "betaCoeff".  The
%     number of rows in attr, compAttr, and betaCoeff must match.  The first row
%     corresponds to the price. As a shorthand, the first element of "attr"
%     is cost rather than price.  The return variable "profit" has the same
%     units as attr(1).

%  Algorithm: first, the price that maximizes the profit is computed. The
%  solution to this optimization problem can be expressed analytically
%  using the inverse Lambert W-function, which in turn can be computed
%  very efficiently with only 3 to 4 Newton-Raphson iterations.


% check the validity of the inputs
if nargin<4
    market = 1e5; % if not provided, use the default value of 100,000
    if nargin<3
        % If not provided, use as a default model the model for a battery
        % electric vehicle, with price in [$], range in [km],
        % 0-100km/h-time in [s],and maximum velocity in [km/h].  The
        % coefficients for this logit model were obtained from student
        % surveys performed in the past.
        betaCoeff = [-1; 40; -2800; 120]/1000;
        if nargin<2
            % if not provided, use the default competition attributes:
            % price of $40k, range of 200km, 0-100km/h-time of 6s, and
            % maximum velocity of 220 km/h. 
            compAttr = [55000;200;6;220];  
        end
    end
end

[n,m] = size(betaCoeff);
if m>1
    error('betaCoeff must be column vector');
end

[na,ma] = size(attr);
if (n~=na || m~=ma)
    error('attr must be the same size as betaCoeff');
end

if (size(compAttr,1) ~= n)
    error('number of columns in attr must match compAttr');
end

% compute the offset resulting from the competition
expVcomp = sum(exp(betaCoeff'*compAttr));
invLambert = exp(betaCoeff'*attr - 1)/expVcomp;

% Iteratively solve for the inverse Lambert function for the input value
% invLambert
guess = 1;
prevGuess = 0;
while abs((guess-prevGuess)/guess) > 1e-12
    prevGuess = guess;
    guess = guess - (guess*exp(guess)-invLambert)...
          /(exp(guess)*(guess+1)-(guess+2)*(guess*exp(guess)-invLambert)/(2*guess+2));
end

% optimal price
price = -(guess-betaCoeff(1)*attr(1)+1)/betaCoeff(1);

% given the price, compute the corresponding demand and profit
demand = market*exp(betaCoeff'*[price;attr(2:end)])...
       /(exp(betaCoeff'*[price;attr(2:end)])+expVcomp);
profit = (price - attr(1))*demand;
end

