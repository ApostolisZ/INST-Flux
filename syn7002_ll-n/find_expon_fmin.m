<<<<<<< HEAD
clear
load('labeling_data');
mids = ldata.mids;
mids_std = ldata.mids_std;
t = ldata.t';
numfigs = length(mids);
dimfigs= ceil(sqrt(numfigs));

%find best sum of exponentials that describes the mids; MIDS =
%sum(a_i*exp(b_i*t))

%set the number of exponentials
numvars = 10;
tic;
% maxiter = [1e3 1e4 1e5 1e6];
% for id = 1:4
figure
hold on
for imetab = 1:length(mids)
    mid_full = mids{imetab};
    subplot(dimfigs,dimfigs,imetab)
    hold on
    for istate = 1:size(mid_full,1)
        mids_1 = mids{imetab}(istate,:);
        mids_std_1 = mids_std{imetab}(istate,:);
        
        % options = optimset('fminsearch');
        % optnew = optimset(options,'MaxFunEvals',1e6,'MaxIter',maxiter(id));
        options = optimset('GradObj','on');
        % options = optimset(options,'GradObj','true');
        %get random variables
        for j = 1:1000
            ab0{j} = (rand(2*numvars,1)-rand(2*numvars,1))*0.1;
        end
        
        % parfor (i = 1:12,4)
        parfor (i = 1:1000,8)
            [abfin{i},fval{i},exitflag,output] = fminunc(@(ab) exponential_sim_error(ab,t,mids_1,mids_std_1),ab0{i},options);
        end
        [minfval,idx] = min(cell2mat(fval));
        minfvals{imetab}{istate} = minfval;
        mid_params{imetab}{istate} = abfin{idx};
        a = mid_params{imetab}{istate}(1:numvars);
        b = mid_params{imetab}{istate}(numvars+1:end);
        b = reshape(b,1,[]);
        M = t*b;
        EM = exp(M);
        a = reshape(a,[],1);
        sim_mid = EM*a;
        scatter(t,mids_1)
        plot(sim_mid)
    end
end
toc
% results{id}{1} = minfvals;
% results{id}{2} = mid_params;
% results{id}{3} = toc;
% end



=======
clear
load('labeling_data');
mids = ldata.mids;
mids_std = ldata.mids_std;
t = ldata.t';
numfigs = length(mids);
dimfigs= ceil(sqrt(numfigs));

%find best sum of exponentials that describes the mids; MIDS =
%sum(a_i*exp(b_i*t))

%set the number of exponentials
numvars = 10;
tic;
% maxiter = [1e3 1e4 1e5 1e6];
% for id = 1:4
figure
hold on
for imetab = 1:length(mids)
    mid_full = mids{imetab};
    subplot(dimfigs,dimfigs,imetab)
    hold on
    for istate = 1:size(mid_full,1)
        mids_1 = mids{imetab}(istate,:);
        mids_std_1 = mids_std{imetab}(istate,:);
        
        % options = optimset('fminsearch');
        % optnew = optimset(options,'MaxFunEvals',1e6,'MaxIter',maxiter(id));
        options = optimset('GradObj','on');
        % options = optimset(options,'GradObj','true');
        %get random variables
        for j = 1:1000
            ab0{j} = (rand(2*numvars,1)-rand(2*numvars,1))*0.1;
        end
        
        % parfor (i = 1:12,4)
        parfor (i = 1:1000,8)
            [abfin{i},fval{i},exitflag,output] = fminunc(@(ab) exponential_sim_error(ab,t,mids_1,mids_std_1),ab0{i},options);
        end
        [minfval,idx] = min(cell2mat(fval));
        minfvals{imetab}{istate} = minfval;
        mid_params{imetab}{istate} = abfin{idx};
        a = mid_params{imetab}{istate}(1:numvars);
        b = mid_params{imetab}{istate}(numvars+1:end);
        b = reshape(b,1,[]);
        M = t*b;
        EM = exp(M);
        a = reshape(a,[],1);
        sim_mid = EM*a;
        scatter(t,mids_1)
        plot(sim_mid)
    end
end
toc
% results{id}{1} = minfvals;
% results{id}{2} = mid_params;
% results{id}{3} = toc;
% end



>>>>>>> 2024247f7ad30fed5e7c393c7eddc78b1a8690d7
