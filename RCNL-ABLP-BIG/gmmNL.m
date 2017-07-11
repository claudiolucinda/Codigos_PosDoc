function [f,df] = gmmNL(theta2,Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Esta função calcula o valor da função objetivo GMM
% Adaptada e otimizada por Cláudio R. Lucinda
% FGV-EESP e EAESP
% 2008/09
% Dois Outputs:
% f - Valor da função objetivo
% df - Gradiente da Função - Para otimizadores que precisam
%
% Inputs:
% theta2 - Vetor de Parametros (o theta2w "amassado")
% Data - lots of other stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames=fieldnames(Data);
for i=1:length(fnames)
    eval([fnames{i} '=Data.' fnames{i} ';']);
end

% Carregando
delta = deltaNL_PAR(theta2,Data);

% Regular Version
% delta = deltaNL_PAR(theta2,Data);


% O que fazer quando o algoritmo sai fora
if max(isnan(delta)|isinf(delta))  == 1
    f = 1e+10;
    gmmresid=(1e+10)*ones(size(delta));
    save([data_dir 'gmmresid.mat'],'gmmresid')
    disp(['GMM:  ' num2str(f)])
    if nargout>1
        df=((1e+10)*ones(size(theta2)));
            disp('Gradient                Theta                  ');
        disp('-----------------------------------------------');
        disp([df theta2])
    else
        disp('   Theta              ');
        disp('----------------------');
        disp(theta2)
    end

    
else
    warning('off','MATLAB:nearlySingularMatrix');
    
%    W = (IV'*IV) \ eye(size(IV,2));
    theta1=(x1'*IV * W * IV'*x1)\(x1'*IV * W * IV'*delta);
    gmmresid = delta - x1*theta1;
    f=(gmmresid'*IV)*W*(gmmresid'*IV)';
    save([data_dir 'gmmresid.mat'],'gmmresid')
    disp(['GMM:  ' num2str(f)])
    if nargout>1
        % % Parte do gradiente
        if max(isnan(delta)) ~= 1
            load([data_dir 'mvaloldNL.mat'],'delta0NL');
            temp = jacobNL(delta0NL,theta2,Data)';
            format long e
            df = 2*temp*IV*W*IV'*gmmresid;
        end
        disp('Gradient                Theta                  ');
        disp('-----------------------------------------------');
        disp([df theta2])
    else
        disp('   Theta              ');
        disp('----------------------');
        disp(theta2)
    end

    
end


% 
% if nargout>1
%     % % Parte do gradiente
%     if max(isnan(delta)) ~= 1
%         load([data_dir 'mvaloldNL.mat'],'delta0NL');
%         temp = jacobNL(delta0NL,theta2,Data)';
%         format long e
%         df = 2*temp*IV*W*IV'*gmmresid;
%     end
%     disp('Gradient                Theta                  ');
%     disp('-----------------------------------------------');
%     disp([df theta2])
% else
%     disp('   Theta              ');
%     disp('----------------------');
%     disp(theta2)
% end

    
