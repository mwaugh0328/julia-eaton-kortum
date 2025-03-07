function [m, rec_low_price, rec_cost, issecond, istraded] = sim_trade_pattern_BEJKe(S,tau,theta,sigma,code)%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to simmulate a pattern of trade and then generate a trade
% share matrix and a random sample of final goods prices given bertrand
% pricing from BEJK (2003). 
%
% Note the theta here is the regular one, not AL value 1/theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for goods and countries and the sample size for prices

Ngoods = 100000; % Adjust this number if it is running really slow (not too low though).
Ncntry = length(S);

% Parameters for technologies
eta = sigma;
inveta = 1./(1- eta);
invNgoods = 1./Ngoods;

markup = eta./(eta-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw productivities and and compute unit costs to produce each good in
% each country
rng(03281978+code)

p1const = zeros(Ngoods,Ncntry);
p2const = zeros(Ngoods,Ncntry); 

for j = 1: Ncntry
    
    u1 = rand(Ngoods,1);
   
    
    z_one = (log(u1)./(-S(j))).^(-1./theta); 
    % First draw, max productivity, this comes from 1-exp(-Sz_one^-theta)   
        
    p1const(:,j) = z_one.^-1; % Invert to convert to unit costs
        
    % Second draw, second best productivity, this comes from
    % 1-exp(-Sz_two^-theta +Sz_one^-theta) 
    
end

for j = 1: Ncntry

 u2 = rand(Ngoods,1);
 
 z_two = (log(u2)./(-S(j)) + (p1const(:,j).^-1).^(-theta)).^(-1./theta);
 
 p2const(:,j) = z_two.^-1;
 
end


% clear u1 z_one u2 z_two % clear large variables to clear up memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop to calculate the low price and country suppliers
%

m = zeros(Ncntry,Ncntry);
sum_price = zeros(Ncntry,1);
rec_low_price = zeros(Ngoods,Ncntry); 
rec_cost = zeros(Ngoods,Ncntry); 
issecond = zeros(Ngoods,Ncntry);
istraded = zeros(Ngoods,Ncntry);


for gd = 1:Ngoods  % This is the good
    
    for im = 1:Ncntry  % This is the country importing the good
         
         % The issue here relative to the EK code is we need to keep track
         % of the second low cost supplier. The sort command in matlab is
         % very slow and time consumeing so here is the strategy:
         
         % First compute the first two entries of potential prices and
         % order them, lowest and second lowest
        
         cif_price = zeros(Ncntry,1);
        
         cif_price(1) = tau(1,im).*p1const(gd,1);
         cif_price(2) = tau(2,im).*p1const(gd,2);
        
         if cif_price(1)< cif_price(2)
            low_cost = cif_price(1);
            low2_cost = cif_price(2);
         else
            low_cost = cif_price(2);
            low2_cost = cif_price(1);
         end
        
            for ex = 3:Ncntry   % This is the country (potentially) exporting the good

                cif_price(ex) = tau(ex,im).*p1const(gd,ex);
                               
                % Now we need to figuer out the potentially new second low
                % cost supplier, so compare the largest guy from above to
                % the second low cost supplier and the minnimum is the new
                % second low cost supplier.
                                
                low2_cost = min(max(cif_price(ex),low_cost),low2_cost);
                
                % Now check if ex is the low cost supplier and record it
                
                low_cost = min(cif_price(ex),low_cost);
            end
            price_charged = 0;

            % First check if home guy is low cost
            if (low_cost == cif_price(im)) 
                istraded(gd,im) = 0;
                
                price_charged = min( min( p2const(gd,im) , low2_cost ) , markup.*low_cost);
                % There is no tau on p2const b.c. this is the home/home
                % situation 
                
                m(im,im) = m(im,im) + price_charged.^(1-eta);
            else
                istraded(gd,im) = 1;
            % If not, then keep looking across all countries    
            for ex = 1:Ncntry	
                
                if (low_cost == cif_price(ex))
                    
                    % Like EK, find the low cost guy, record trade share.
                    
                    price_charged = min( min( tau(ex,im).*p2const(gd,ex) , low2_cost ) , markup.*low_cost);
                                                    
                    % This is key here:
                    % First take the mimium cost between the second best domestic 
                    % competitor or the exporter and the second lowest cost competitors, this is 
                    % globally the best guy. (Might want to walk through this logic as to why)

                    % Then the low cost guy will not want to charge more than a
                    % markup (eta/(eta-1)over marginal cost, so it is either that
                    % or if lower,i.e. the guy limit prices against the competitor.
            
                    m(ex,im) = m(ex,im) + price_charged.^(1-eta);
                    break
                    
                    % Note that the price is not low_cost, but some markup
                    % over marginal cost.
                    % Note previous code exploted eta = 2, so I had 1/price
                    % rather then this. I think this is correct, see b.52
                    % in SW (2011).                    
                end
            end

            end
        
        % Now record the low cost guy, but price charge again.

        sum_price(im) = sum_price(im) + price_charged.^(1-eta); 
        % Note previous code exploted eta = 2, so I had 1/price
        % rather then this. I think this is correct
        
        rec_low_price(gd,im) = price_charged;
        rec_cost(gd,im) = low_cost;
        
        if price_charged == low2_cost
            issecond(gd,im) = 1;
        else
            issecond(gd,im) = 0;
        end
                
        % rec_low_price(gd,im) = low_cost; Use this line, it will then look
        % exactly like EK when the method of moment estimation is applied.
        

%         rec_markup(gd,im) = price_charged./low_cost;
        % As a test, record markups. This should be pareto on the interval
        % 1 to monopoly markup m

    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop to calculate aggregate price index and the trade shares.

g = zeros(Ncntry,1);

for im = 1:Ncntry										

	g(im) = (sum_price(im)*invNgoods).^(inveta);
															
    for ex = 1:Ncntry
    
        m(ex,im) = invNgoods.*m(ex,im)./g(im).^(1- eta);    
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a price matrix simmilar to that in the data

% rand('twister',02071983+code)
% 
% keep = randi(Ngoods,sample,1);
% include = sample./Ngoods; % Vary the sample to see how it affects the estimate  of theta  
% keep = u < include;

% Take the prices from each 

% final_cost = rec_low_cost;
% final_markup = rec_markup;




