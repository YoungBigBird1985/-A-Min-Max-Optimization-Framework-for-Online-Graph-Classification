function run_ridge_min_max_release()

addpath('D:\SSLGC algorithm\data');

dataname = ('run_DAGGER_ridge_min_max');
coraLabelNum = 7;
Pubmed_DiabetesLabelNum = 3;
IMDB = 4;
coauthor = 4;

resultFileName = [dataname '.log'];
fid = fopen(resultFileName,'a');
result_description = ('sample No. error rate query No. computational time' );

load('Coauthor.mat');

LableM_R = gnd; % 2485*1
M_lr = M'; % M: 2485*100
Dimention = size(LableM_R,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LowRank = 101;
NumberShuffl = 21;
c = 0.01;
history_record_err_all = []; history_record_query_all = []; history_record_time_all = [];
for c_index = -3:1
        mu = 10^(c_index); mu
        errorrate_all = []; querynumber_all = [];
        errorrate_class = []; querynumber_class = []; updatenumber_class = []; time_class = [];
        history_record_err_collect = []; history_record_query_collect = []; history_record_time_collect = [];
for labelindex = 1:coauthor
    LableM = zeros(Dimention,1);
    for i = 1:Dimention
        if(LableM_R(i,1) == labelindex)
            LableM(i,1) = 1;
        else
            LableM(i,1) = -1;
        end
    end

    errorrate_avg = []; querynumber_avg = []; updatenumber_avg = []; time_avg = [];
    for j=1:NumberShuffl
        
        index = ID_ALL(j,:);
        M_lr_shufl = M_lr(:, index);
        LableM_shufl = LableM(index,:);
        
        errorcount = 0; error_rate = 0;
        querynumber = 0; updatenumber = 0;   

        A = mu.*eye(LowRank-1); b = zeros(LowRank-1,1); w = zeros(LowRank-1,1);
        inv_A = mu^(-1).*eye(LowRank-1);
        
        Q_t = 0; Z_t = 0; 
        history_record_err = []; history_record_query = []; history_record_time = [];
        
        tic; t1 = clock;
        for t = 1:Dimention
            m_t = M_lr_shufl(:,t);
            inv_A_add = inv_A - (inv_A*m_t*m_t'*inv_A)./(0.5+m_t'*inv_A*m_t);
            P_t = m_t'*inv_A_add*b;
            Y_P = sign(P_t); % sign(w(T)m)
            Y_A = LableM_shufl(t,1);
            r = m_t'*inv_A*m_t;
            F = P_t*P_t + 2*abs(P_t) - r/(1+2*r);
            if  F >= 0
                %draw a nernoulli random variable Q with probability 
                Bernoulli_P = 2*c/(2*c + F);
                Q_t = binornd(1,Bernoulli_P);
                 if Q_t == 1  
                      querynumber = querynumber + 1;
                      if Y_P ~= Y_A
                          Z_t = 1;
                      end
                 end
            else
                Q_t = 1; Z_t = 1;
                querynumber = querynumber + 1;
            end
            
            if Y_P ~= Y_A
                errorcount = errorcount + 1;
            end
            if Q_t == 1 && Z_t == 1
                 updatenumber = updatenumber + 1;
                 A = A + 2*m_t*m_t';
                 b = b + Y_A.*m_t;
                 inv_A = inv_A - (inv_A*m_t*m_t'*inv_A)./(0.5+m_t'*inv_A*m_t); % Sherman-Morrisan Identity
                 %w = inv_A*b;
            end
            Q_t = 0; Z_t = 0;
            
            error_rate = errorcount/t;
            if mod(t, 50) == 1
                history_record_time = [history_record_time, etime(clock,t1)];
                history_record_query = [history_record_query, querynumber];;
                history_record_err = [history_record_err, error_rate];
            end
        end % t dimention index
        errorrate_avg = [errorrate_avg, error_rate];        errorrate_all = [errorrate_all, error_rate];
        querynumber_avg = [querynumber_avg, querynumber];   querynumber_all = [querynumber_all, querynumber];
        updatenumber_avg = [updatenumber_avg, updatenumber];
        time_avg = [time_avg, etime(clock,t1)];
        history_record_time_collect = [history_record_time_collect; history_record_time]; history_record_time = [];
        history_record_query_collect = [history_record_query_collect; history_record_query]; history_record_query = [];
        history_record_err_collect = [history_record_err_collect; history_record_err]; history_record_err = [];
    end
    %average_result = ['average error rate: ', num2str( mean(errorrate_avg) ), ' average query number: ', num2str(mean(querynumber_avg)), ' average time cost: ', num2str(mean(time_avg))];
    %fprintf(fid, '%s\n', average_result);

    labelindex
    errorrate_class = [errorrate_class, mean(errorrate_avg)];
    querynumber_class = [querynumber_class, mean(querynumber_avg)];
    updatenumber_class = [updatenumber_class, mean(updatenumber_avg)];
    time_class = [time_class, mean(time_avg)];
% end  % rank for
end  % label index for
fprintf('error rate %g, %g\n', mean(errorrate_class), std(errorrate_all, 1));
fprintf('query number %g, %g\n', mean(querynumber_class), std(querynumber_all, 1) );
%mean(updatenumber_class)
%mean(time_class)
history_record_time_all = [history_record_time_all; mean(history_record_time_collect)];
history_record_query_all = [history_record_query_all; mean(history_record_query_collect)];
history_record_err_all = [history_record_err_all; mean(history_record_err_collect)];
end % mu index for

x = 1:50:Dimention;
plot(x, history_record_err_all);
for i=1:size(history_record_err_all,1)
    fprintf(fid, '%g\t', history_record_err_all(i,:));
    fprintf(fid, '\n');
end

x = 1:50:Dimention;
plot(x, history_record_err_all); hold on
plot(x, history_record_query_all); hold on
plot(x, history_record_time_all);
fprintf(fid, 'error history record\n');
for i=1:size(history_record_err_all,1)
    fprintf(fid, '%g\t', history_record_err_all(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, 'query history record\n');
for i=1:size(history_record_query_all,1)
    fprintf(fid, '%g\t', history_record_query_all(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, 'time history record\n');
for i=1:size(history_record_time_all,1)
    fprintf(fid, '%g\t', history_record_time_all(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

