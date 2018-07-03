clear
load('labeling_data');
syms ts
mids = ldata.mids;
mids_std = ldata.mids_std;
t = ldata.t';
numfigs = length(mids);
dimfigs= ceil(sqrt(numfigs));
%find best sum of exponentials that describes the mids; MIDS =
%sum(a_i*exp(b_i*t))

%set the number of exponentials
numvars =2;
tic;
maxiter = [1e3 1e4 1e5 1e6];
% for id = 1:4
figure

for imetab = 1:length(mids)
    mid_full = mids{imetab};
    subplot(dimfigs,dimfigs-1,imetab)
%     figure
    hold on
    for istate = 1:size(mid_full,1)
        mids_1 = mids{imetab}(istate,:);
        mids_std_1 = mids_std{imetab}(istate,:);
        
         options = optimset('Display','off');
        % optnew = optimset(options,'MaxFunEvals',1e6,'MaxIter',maxiter(id));
        
        %get random variables
        for j = 1:1000
            ab0{j} = (rand(2*numvars,1)-rand(2*numvars,1))*10;
        end
        
        parfor (i = 1:1000,8)
            [abfin{i},fval{i},exitflag,output] = fminsearch(@(ab) exponential_sim_error(ab,t,mids_1,mids_std_1),ab0{i},options);
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
        sim_mid_f = sum(a'.*exp(b*ts));
        avg = sum(mids_1)/length(mids_1);
        sstot = sum((mids_1-avg).^2);
        ssres = sum((mids_1-sim_mid').^2);
        rsq{imetab}(istate) = 1-ssres/sstot;
        fprintf([num2str(rsq{imetab}(istate)) '  '])
%         plot(t,sim_mid)
        fplot(sim_mid_f,[t(1),t(end)]);
%         legend([repmat(['M + '], istate,1) num2str((0:istate-1)')])
        scatter(t,mids_1)

    end
    [worst_rsq(imetab),worst_idx(imetab)] = min(rsq{imetab});
%      legend([repmat(['M + '], istate,1) num2str((0:istate-1)')])
    fprintf('\n');
end
toc
foo = 1;

% results{id}{1} = minfvals;
% results{id}{2} = mid_params;
% results{id}{3} = toc;
% end



