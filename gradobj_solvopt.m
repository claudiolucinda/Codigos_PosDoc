function df = gradobj_solvopt(theta2)

% This function computes the gradient of the objective function
% Based on code written by Aviv Nevo, May 1998
% Adapted by Matthijs Wildenbeest, April 2010
% Adapted by Claudio Lucinda, Sept 2013

global args

fnames=fieldnames(args);
for i=1:length(fnames)
    eval([fnames{i} '=args.' fnames{i} ';']);
end

%args.cmg=mc.*(1-not_ok)-ones(size(mc)).*not_ok;
miolo=pinv(IV);
gmmresid = load([data_dir '\gmmresid.mat']);    
%mval1 = load([data_dir '\mvalold.mat']);
temp = jacob_solvopt(theta2)';
format long e

df = 2*temp*IV*miolo*gmmresid.gmmresid;