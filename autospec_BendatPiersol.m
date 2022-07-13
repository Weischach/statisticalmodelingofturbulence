close all; clear all; clc;





%% There are 101 files with wind data from multiple anemometers
%% filenames.mat contains the name of these 101 files
%% mnum=9 will load the file 20141031-Modint-B02.mat
load filenames.mat

for mnum = 9
    
    %% load the experimental data containing wind measured in anemometers  
    m = fname{mnum};
    load(m);
    
    
    %% filter the wind data to remove duration with high wind angle at the beginning and end
    ly = length(yaw);
    sy = floor(ly/40);
    while abs(mean(yaw(1:sy)))>5
        k = 1;
        while abs(yaw(k))>5 && k<ly
            k = k+1;
        end
        fprintf('k = %d', k)
        yaw(1:k) = [];
        w(1:k,:) = [];
        v(1:k,:) = [];
        ly = length(yaw);
        sy = floor(ly/40);
    end
    
    if ly>0
        while abs(mean(yaw(end-sy:end)))>5
            k = 1;
            while abs(yaw(end+1-k))>5 && k<ly
                k = k+1;
            end
            fprintf('k = %d \n', k)
            yaw(end+1-k:end) = [];
            w(end+1-k:end,:) = [];
            v(end+1-k:end,:) = [];
            ly = length(yaw);
            sy = floor(ly/40);
        end
    end
    fprintf('yaw is reduced to length %d for file number %d \n', ly, mnum');
    
    
    %% select 1 for vertical data, and 2 for lateral data. and select column number to select anemometer port.
    for kk1 = 1
        if kk1 == 1
            data = w;
            m1a = strcat(m,'dataw');
        end
        if kk1 == 2
            data = v;
            m1a = strcat(m,'datav');
        end
    
    fprintf('Working on file number %d \n', mnum');
    
    for columnu = 1
        m1b = strcat(m1a,num2str(columnu))
        m1c = strcat(m1b,'.mat')
        %% load m1c if all the datafiles are saved and start the code from there.
        % calculating turbulence and standard deviation from the wind data
        time = time(1:ly);
        dataw = data(:,columnu);
        lentemp = floor(length(dataw)/1);
        dataw = dataw(1:lentemp);
        lenw = length(dataw);
        W = mean(dataw);
        turbw = dataw - W;
        variancew = mean(turbw.^2);
        sigmaw = sqrt(variancew);
        deltat = 0.1
        
        
        %% certain anemometers have ot collected headwind, and hence all the data gets filtered.
        % finding the autospectrum of the wind data
        if mean(abs(turbw))>0
            psdpts = 2^nextpow2(length(turbw)/16);
            psdlim = 5;
            [pww1,fww1] = PSDofRAW(turbw, psdpts, deltat, psdlim);
            pww1 = smooth(pww1);
            fww = [0 logspace(-3, 0.695, 499)];
            pww = interp1(fww1,pww1,fww);
            T_wchar = pww(2)/(4*variancew);
            ft = trapz(fww,pww);
            lt = abs(ft-variancew)*100/variancew;
            fprintf('The difference between std. devn and periodogram integral is %d \n', lt)
            save(m1b, 'turbw', 'time', 'sigmaw', 'T_wchar', 'data');
            
            
            coeffs = 3;
            alpha_coeff = [0.373 0.1996 0.122 0.085 0.064  0.050 0.0409];
            A_coeff = [0.186 0.373 0.558 0.745 0.931 1.117 1.30];
            alpha_coeff = alpha_coeff(1:coeffs);
            A_coeff = A_coeff(1:coeffs);
            
            
            Betas0 = zeros(1,coeffs);
            Betas0(1) = 1;
            Betas = Betas0;
            Aeq = [ones(1,coeffs)];
            beq = [1];
            Aineq = -alpha_coeff;
            bineq = 0.0001;
            tau = [0:0.01:100];     
            
            %% Physical qunatity A can only be found on case by case basis. The optimization can also be run after finding A for the particular case.
            % if A is found, the optimization can be done using wpwwerror.
            % Otherwise wpwwerrornoA can be used.
            
            
            options= optimset('Algorithm','sqp', 'MaxIter', 200, 'TolCon',1e-10,'TolFun',1e-8,'TolX',1e-8,'Diagnostics','off');
            options.MeshTolerance = 1e-5;
            
            
            [Betas, err]=patternsearch(@(Betas)wpwwerrornoA(Betas, pww, fww, alpha_coeff, T_wchar, tau, variancew, coeffs),Betas0,Aineq,bineq,Aeq,beq,-Inf*ones(1,coeffs),Inf*ones(1,coeffs),[],options);
            
            
            %% after the beta coefficients are found, the perturbation series autopsectrum can be derived.
            alpha = dot(Betas,alpha_coeff);
            tuu = (alpha/T_wchar)*tau;
            tuuvon = (alpha_coeff(1)/T_wchar)*tau;
            Ruu = ((2^(2/3))/gamma(1/3))*[(tuu.^(1/3)).*besselk(1/3,tuu) - 0.5*(tuu.^(4/3)).*besselk(2/3,tuu)];
            Ruu(1) = 1;
            Ruuvon = ((2^(2/3))/gamma(1/3))*[(tuuvon.^(1/3)).*besselk(1/3,tuuvon) - 0.5*(tuuvon.^(4/3)).*besselk(2/3,tuuvon)];
            Ruuvon(1) = 1;
            pwwperturb = psdperturbseries(Betas, fww, Ruu, tau, variancew, coeffs);
            pwwvonKarman = psdperturbseries(Betas0, fww, Ruuvon, tau, variancew, coeffs);
            vonKarmanerror = wpwwerrornoA(Betas0, pww, fww, alpha_coeff, T_wchar, tau, variancew, coeffs);
            lt2 = 100*abs(variancew-trapz(fww,pwwperturb))/variancew;
            ltf = 2*(lt^2 + lt2^2);
            
            
            %% the standarddeviation is calculated by integrating the autospectrum, and the autospectra that deviates (according to Equation 8 in paper) are neglected
            if (ltf)>900
                m1 = strcat(m1b,'_PSD512smooth1_variiance_2')
            elseif min(pwwperturb) < 0
                m1 = strcat(m1b,'_PSD512smooth1_neg_2')
            else
                m1 = strcat(m1b,'PSD512smooth1_2')            
                save(m1, 'Betas', 'variancew', 'T_wchar', 'deltat', 'pww', 'fww', 'pwwvonKarman', 'pwwperturb');              
                
                close all;
                fig1 = figure;
                loglog(fww,pww, ':k', fww,pwwperturb, '-k', fww,pwwvonKarman, '--k')
                hold on
                loglog(fww,pww, ':k', 'LineWidth', 0.5)
                loglog(fww,pwwperturb, '-k', 'LineWidth', 3)
                loglog(fww,pwwvonKarman, '--k', 'LineWidth', 1)
                legend('Data', 'Perturbation series', 'von Karman')
                ylabel('Autospectrum S_w_w(f)')
                xlabel('Frequency (f)')
                grid on
                savefig(m1)
                saveas(fig1,m1,'jpeg')
            end
        end
    end
    end
end
