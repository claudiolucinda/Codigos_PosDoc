% File jointly generates income/age/gender/HHcomposition/HHsize draws for whole country Sweden
% Cristian Huse -- July 2011
% --------------------------------------------------------------------------------------------------------------------------------

clear; clc;

% settings
% ----------
factor = 1;             %   factor by which population is divided in order to generate draws
ns = 1000;              %  number of simulations/year
descriptives = 0;   %   =1 for some basic descriptives

% set dir, load file, read, define basic vars.
% --------------------------------------------------
%cd C:/Atual/Green_Cars_Sweden/Data/Flex_fuel_sweden;

%cd 'C:\Documents and Settings\Claudio Lucinda\Meus documentos\Papers\Christian, Salvo\Flex_fuel_sweden';

cd 'C:\Users\claudiolucinda\Documents\Flex_fuel_sweden';

% Distribuição de Renda
[num,txt,raw] = xlsread('draws_jul2011.xls');
size(num); size(txt); size(raw);

% Selecionando os valores da renda
num = num(2:size(num,1),:);
% selecionando os mínimos de cada bracket
mini      =   num(:,1);
% selecionando os máximos de cada bracket
maxi     =   num(:,2);
% selecionando os dados de renda
draw0   =   num(:,3:9);
% selecionando o sexo
gender  =   num(:,10);% gender = 1 if male; 0 otherwise
% selecionando as faixas de idade
age       =   num(:,11);

age_brakts_low=[20;25;30;35;40;45;50;55;60;65;70;75;80;85];
age_brakts_hi=[25;30;35;40;45;50;55;60;65;70;75;80;85;100];
% name dimensions for use in loop
nrows =size(draw0,1);
ncols  =size(draw0,2);

% load second data file
% Age, Household composition & Household size
[num2,txt2,raw2] = xlsread('draws2_jul2011.xls');

% selecionando os valores
num2 = num2(2:size(num2,1),:);
draw1   =   num2(:,4:10);
HHcomp  =   num2(:,11);
HHcompnum   = num2(:,12);
HHsize  =   num2(:,13);

% name dimensions for use in loop2
nrows2 =size(draw1,1);
ncols2  =size(draw1,2);

% initialize vars just in case

%tmp1=[]; tmp2=[]; tmp3=[]; tmp4=[]; tmp5=[]; inc_tmp=[]; gender_tmp=[]; age_tmp=[];
inc_accum03 =[]; age_accum03=[]; gender_accum03=[]; inc_accum04 =[]; age_accum04=[]; gender_accum04=[];
inc_accum05 =[]; age_accum05=[]; gender_accum05=[]; inc_accum06 =[]; age_accum06=[]; gender_accum06=[];
inc_accum07 =[]; age_accum07=[]; gender_accum07=[]; inc_accum08 =[]; age_accum08=[]; gender_accum08=[];
inc_accum09 =[]; age_accum09=[]; gender_accum09=[]; 

% big loop
% ----------
for t = 1:ncols;    % loop across columns i.e.years - Banco dados renda

    tmp1 = draw0(:,t);  % data for year t
            
    if factor > 0 ;         % see settings
        tmp3 = ceil(tmp1./factor) ;
    elseif factor ==0;
        tmp2 =min(tmp1);
        tmp3 = ceil(tmp1./tmp2);
    end
    
    for i = 1:nrows;    % loop across rows i.e. income-age-gender combinatons
        tmp4 = mini(i,:);   % for each income bracket
        tmp5 = maxi(i,:);
        mina=age(i,:);
        if mina~=85
            maxa=mina+5;
        else
            maxa=100;
        end
 
        if tmp3(i,:)>0;
            tmp6 = tmp3(i,:);
        
            inc_tmp = random('unif', tmp4, tmp5,tmp6,1);    % generate income - Uniform Draws from each bracket
        
            gender_tmp = gender(i,:)*ones(tmp6,1);              % retrieve gender & age
            age_tmp = random('unif',mina,maxa,tmp6,1);
        
            eval(['inc_accum0' num2str(t+2) '=[inc_accum0' num2str(t+2) ';inc_tmp];']);
            eval(['gender_accum0' num2str(t+2) '=[gender_accum0' num2str(t+2) ';gender_tmp];']);
            eval(['age_accum0' num2str(t+2) '=[age_accum0' num2str(t+2) ';age_tmp];']);
        end 
        
    end

end

% Testando se os draws de aleatório são melhores ou piores do que os Halton
% draws

% select = random('unid', size(inc_accum03,1),ns,1);
% tmp=tiedrank(age_accum03)/length(age_accum03);
% tmp2=tiedrank(inc_accum03)/length(inc_accum03);
% ttmp=round(1e4*[tmp tmp2])/1e4;
% p = haltonset(2,'Skip',1e3,'Leap',1e2);
% p = scramble(p,'RR2');
% select2=net(p,10*ns);
% % arredondando em 10e4
% select2=round(1e4*select2)/1e4;
% a=ismember(select2,ttmp,'rows');
% test1=age_accum03(select,:);
% test2=age_accum03(b,:);


%initialize more vars.
HHcompnum03=[];  HHcompnum04=[];  HHcompnum05=[];  HHcompnum06=[];  HHcompnum07=[];  HHcompnum08=[];  HHcompnum09=[];
HHsize03=[];  HHsize04=[];  HHsize05=[];  HHsize06=[];  HHsize07=[];  HHsize08=[];  HHsize09=[];

for t = 1:ncols2;    % loop across columns i.e.years

    tmp1 = draw1(:,t);  % data for year t
            
    if factor > 0 ;         % see settings
        tmp3 = ceil(tmp1./factor) ;
    elseif factor ==0;
        tmp2 = ( min(tmp1) );
        tmp3 = ceil(tmp1./tmp2);
    end;
    
    for i = 1:nrows2;    % loop across rows i.e. income-age-gender combinatons
 
        if tmp3(i,:)>0;
            tmp6 = tmp3(i,:);

            HHcompnum_tmp   = HHcompnum(i,:)*ones(tmp6,1);
            HHsize_tmp  =   HHsize(i,:)*ones(tmp6,1);

            eval(['HHcompnum0' num2str(t+2) '=[HHcompnum0' num2str(t+2) ';HHcompnum_tmp];']);
            eval(['HHsize0' num2str(t+2) '=[HHsize0' num2str(t+2) ';HHsize_tmp];']);
        end
        
    end

end





% extract ns draws
% --------------------

tresh_inc=500;
tresh_age=65;

for t=1:ncols2
    
    g1=eval(['find(inc_accum0' num2str(t+2) '<tresh_inc);']);
    g2=eval(['find(age_accum0' num2str(t+2) '>tresh_age);']);
    g3=union(g1,g2);
    eval(['inc_accum0' num2str(t+2) '(g3,:)=[];']);
    eval(['gender_accum0' num2str(t+2) '(g3,:)=[];']);
    eval(['age_accum0' num2str(t+2) '(g3,:)=[];']);
    
    
    tempvar=eval(['inc_accum0' num2str(t+2) ';']);
    select=random('unid',size(tempvar,1),ns,1);
    eval(['draws0' num2str(t+2) '=[inc_accum0' num2str(t+2) '(select,:) gender_accum0' num2str(t+2) '(select,:) age_accum0' num2str(t+2) '(select,:)];']);

    select2 = random('unid', size(HHcompnum03,1),ns,1);
    eval(['draws0' num2str(t+2) '=[draws0' num2str(t+2) ' HHcompnum0' num2str(t+2) '(select2,:) HHsize0' num2str(t+2) '(select2,:)];']);
    

end
% checking some descriptive stats
% ----------------------------------------
if descriptives ==1;

tmp = [mean(draws03); std(draws03); range(draws03); quantile(draws03,.05); quantile(draws03,.25); quantile(draws03,.5); quantile(draws03,.75); quantile(draws03,.95)];
corr(draws03);

tmp = [mean(draws04); std(draws04); range(draws04); quantile(draws04,.05); quantile(draws04,.25); quantile(draws04,.5); quantile(draws04,.75); quantile(draws04,.95)];
corr(draws04);

tmp = [mean(draws05); std(draws05); range(draws05); quantile(draws05,.05); quantile(draws05,.25); quantile(draws05,.5); quantile(draws05,.75); quantile(draws05,.95)];
corr(draws05);

tmp = [mean(draws06); std(draws06); range(draws06); quantile(draws06,.05); quantile(draws06,.25); quantile(draws06,.5); quantile(draws06,.75); quantile(draws06,.95)];
corr(draws06);

tmp = [mean(draws07); std(draws07); range(draws07); quantile(draws07,.05); quantile(draws07,.25); quantile(draws07,.5); quantile(draws07,.75); quantile(draws07,.95)];
corr(draws07);

tmp = [mean(draws08); std(draws08); range(draws08); quantile(draws08,.05); quantile(draws08,.25); quantile(draws08,.5); quantile(draws08,.75); quantile(draws08,.95)];
corr(draws08);

tmp = [mean(draws09); std(draws09); range(draws09); quantile(draws09,.05); quantile(draws09,.25); quantile(draws09,.5); quantile(draws09,.75); quantile(draws09,.95)];
corr(draws09);

%corr(draws03)
%   income     gender     age         HHcomp HHsize
%    1.0000    0.3032   -0.0011   -0.0733    0.0033
%    0.3032    1.0000   -0.0088   -0.0186   -0.0000
%   -0.0011   -0.0088    1.0000   -0.0115   -0.0368
%   -0.0733   -0.0186   -0.0115    1.0000    0.6458
%    0.0033   -0.0000   -0.0368    0.6458    1.0000

%corr(draws09)
%1.0000    0.2020    0.0195    0.0606    0.0912
%0.2020    1.0000    0.0700   -0.0141    0.0209
%0.0195    0.0700    1.0000    0.0327    0.0080
%0.0606   -0.0141    0.0327    1.0000    0.8696
%0.0912    0.0209    0.0080    0.8696    1.0000

%Note generally low correlations b/w income and other variables -- only gender seems to play a moderate role, so might be worth drawing jointly (income,gender)
%Note also that this goes against the "household" definition, which might make sense in a country with many singles
%From the second set of draws, note how HHcomp (assigning numbers 1-12 to different household compositions) correlates highly with HHsize
%One idea for HHcomp would be to generate dummies, but not sure how easy to identify
%Applications: while income draws could be used `a la BLP in a term alpha*ln(y_i - p_j) or directly in the price coefficient, HHsize could be used in the size coefficient

% % some histograms
% hist(draws03(:,1))
% hist(draws03(:,3));
% hist(draws03(:,5));
% hist(draws04(:,1));
% hist(draws05(:,1));
% hist(draws06(:,1));
% hist(draws07(:,1));
% hist(draws08(:,1));
% hist(draws09(:,1));
% 
% %some crosstabs
% plotmatrix(draws03(:,[1 2]));figure(gcf);
% %males tend to earn higher salaries
% plotmatrix(draws03(:,[1 3]));figure(gcf);
% %non-monotonic age effect
% plotmatrix(draws03(:,[1 5]));figure(gcf);
% plotmatrix(draws03(:,[2 3]));figure(gcf);
% plotmatrix(draws03(:,[2 5]));figure(gcf);
% plotmatrix(draws03(:,[3 5]));figure(gcf);
% %no clear pattern really
% 
% plotmatrix(draws09(:,[1 2]));figure(gcf);
% plotmatrix(draws09(:,[1 3]));figure(gcf);
% plotmatrix(draws09(:,[1 5]));figure(gcf);
% %broadly consistent w/ the above

end;

% need to deflate
% -------------------

deflat=dlmread('deflat.txt','\t',1,0);

% selecting the deflator from july for each year since 2003

deflindex=deflat((deflat(:,1)>2002 & deflat(:,2)==6),3);

for i=1:7
    j=i+2;
    eval(['d_draws_0' num2str(j) '=[];'])
    for k=1:12
        select = random('unid', size(inc_accum03,1),ns,1);
        eval(['tmp_0' num2str(j) '=inc_accum0' num2str(j) './deflindex(' num2str(i) ');'])
        eval(['tmp_0' num2str(j) '=inc_accum0' num2str(j) '(select,1);'])
        eval(['d_draws_0' num2str(j) '=[d_draws_0' num2str(j) ' tmp_0' num2str(j) '];']);
    end
    
end

for i=1:7
    j=i+2;
    eval(['d_age_0' num2str(j) '=[];'])
    for k=1:12
        tmp=eval(['tiedrank(age_accum0' num2str(j) ')/length(age_accum0' num2str(j) ');']);
        tmp=round(1e4*tmp)/1e4;
        p = haltonset(1,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'RR2');
        select=round(1e4*net(p,ns))/1e4;
        % arredondando em 10e4
        [a,b]=ismember(select,tmp);

        %select = random('unid', size(age_accum03,1),ns,1);
        eval(['tmp_0' num2str(j) '=age_accum0' num2str(j) '(b,1);'])
        eval(['d_age_0' num2str(j) '=[d_age_0' num2str(j) ' tmp_0' num2str(j) '];']);
    end
    
end


d_draws_03=d_draws_03';
d_draws_04=d_draws_04';
d_draws_05=d_draws_05';
d_draws_06=d_draws_06';
d_draws_07=d_draws_07';
d_draws_08=d_draws_08';
d_draws_09=d_draws_09';


d_age_03=d_age_03';
d_age_04=d_age_04';
d_age_05=d_age_05';
d_age_06=d_age_06';
d_age_07=d_age_07';
d_age_08=d_age_08';
d_age_09=d_age_09';

% write output into file
% --------------------------
save('draws_light4.mat','d_draws_03','d_draws_04','d_draws_05','d_draws_06','d_draws_07','d_draws_08','d_draws_09');
save('age_light4.mat','d_age_03','d_age_04','d_age_05','d_age_06','d_age_07','d_age_08','d_age_09');