clear all;
%rng(5);

NITERATIONS=10;

H=50;
W=70;
F=100;
fps=10;

S=2;

Gamma=zeros(NITERATIONS*H,S*W); % Keep track of the evolution of gamma
Qs=zeros(NITERATIONS*W,F,3); %evolution of q

offset=15;
width=20/2;
yfun=@(x)offset+(W-2*offset)*((2*x-F)/F).^2;

%plot(yfun(1:F));

width2=10/2;

frames=.7*ones(H,W,F);
appearance = (rand(F,1)<.6);
for f=1:F,
    if appearance(f), %first appearance
        frames(round(1*end/3):round(2*end/3),round(yfun(f)-5):round(yfun(f)-5+2*width),f)=.1;
        FRAME1=f;
    else
        frames(round(2*end/4):round(3*end/4),round(yfun(f)-5):round(yfun(f)-5+2*width2),f)=.1;
        FRAME2=f;
    end
end

SNRv = sum(frames(:).^2);
noise = 1/10*randn(H,W,F);
frames = frames + noise;
frames(frames<0)=0;
frames(frames>1)=1;
SNRv = 10 * log(SNRv / sum(noise(:).^2));
fprintf('Video SNR: %f dB \n', SNRv);

a=audioinfo('Mseq3D.wav');
Fs=a.SampleRate;
% centrage de l'échantillon sonore pour avoir quelque chose
% d'intéressant à écouter...
start=floor((a.TotalSamples - Fs/fps*F)/2);
stop=floor((a.TotalSamples + Fs/fps*F)/2);
[audioL,Fs] = audioread('UnVoyageCythre.wav', [start stop-1]);
audioL=sum(audioL,2); % on convertit en mono

audioR = audioL;
audio = audioL;
lambda=.9; % coefficient d'atténuation relative
audio(1:numel(audioR),2) = audioR*lambda;

%regrouper l'audio en petit paquets (1 par image)
sam=(Fs/fps);
audio=reshape(audio,[sam,F,2]); % audio(:,frame,channel)

%Décaler les petits paquets selon le modèle...
% tau = alpha l + beta
alpha = .2;
beta = -14;
% delta_samples = Fs*tau

for f=1:F,
    audio(:,f,2)=circshift(audio(:,f,2), round(alpha*yfun(f)+beta));
end

% ajouter du bruit
SNRa = sum(audio(:).^2);
noise = 1/200*randn(sam,F,2);
audio = audio + noise;
SNRa = 10 * log(SNRa / sum(noise(:).^2));
fprintf('Audio SNR: %f dB \n', SNRa);

%plot(audio(:,1));
%hold all
%plot(audio(:,2));

%% EM
% initialisation

% calculer le décalage temporel maximal pour chaque frame. sau=tau*Fs (sau c'est tau exprimé en samples)
sau = zeros(F,1);
for f=1:F,
    [xc,lags]=xcorr(audio(:,f,2),audio(:,f,1));
    [m,i]=max(xc);
    sau(f)=lags(i);
end

%mu_b = .5*ones(H,W);%mean(frames,3); % backgroud is mean of images
mu_b=mean(frames,3);
mu=zeros(H,W,S);
gamma=zeros(H,W,S);
%randF=30;
%FRAMESFIX=[FRAME1, FRAME2];
for s=1:S,
    mu(:,:,s) = frames(:,:,randi(F)); % foreground to each class set to random image
   % mu(:,:,s) = frames(:,:,FRAMESFIX(s));
    gamma(:,:,s) = abs(mu(:,:,s)-mu_b);
    gamma(:,:,s) = gamma(:,:,s)>.5*max(max(gamma(:,:,s)));
%    hpos = mean(gamma(:,:,s),1); 
%    l_back = round(mtimes(1:W,hpos') / sum(hpos)); % translate apparences and masks back
%    mu(:,:,s)=circshift(mu(:,:,s),[0 -l_back]);
%    gamma(:,:,s)=circshift(gamma(:,:,s),[0 -l_back]);
end
%gamma=logical(rand(H,W,S)>.5);
gamma=logical(gamma);


%% Recenter gammas, mu and q
q=ones(W,S,F)/(W*S); % peut importe, on va le changer après
%lq = log(q); % log q(l,s,F)

[ gamma, mu, q ] = recenter( gamma, mu, q, 10 );

O = zeros(F,1); % no initial occlusion
Pi = ones(S,1)/S; % equiprobable apparences

%vsigma = 1/(mean(mean(frames(:,:,1)))); % dummy video inverse covariance
vsigma=10;

asigma=10;

% alpha - beta initial guess:
alpha = .25; % order of magnitude
beta = -W/2*alpha; % assumed to be centered in the center of the screen
A = audio(:,:,1); % true audio channel initial guess
lambda = 1.; % lambda attenuation factor initial guess

%%
for MASTER=1:NITERATIONS,
    
%gamma logging
for s=1:S,
    Gamma((MASTER-1)*H+1:MASTER*H,(s-1)*W+1:s*W)=gamma(:,:,s);
end
lQs=cat(3,squeeze(q(:,1,:)),squeeze(q(:,2,:)),zeros(W,F));
if NITERATIONS>1,
    lQs(lQs==0)=+Inf;
    lQs(lQs==+Inf)=min(lQs(:));
end
lQs=log(lQs);
scaledlQs = (lQs-min(lQs(:))) ./ (max(lQs(:)-min(lQs(:))));
Qs((MASTER-1)*W+1:MASTER*W,1:F,:)=scaledlQs;

    %% E-Step

fprintf('E-step %i\n',MASTER);

lq = zeros(W,S,F);
for l=1:W,
    if mod(l,round(W/5))==0,
         fprintf('translation %i \n',l)
    end
   % q(l,:,:)=1; %l,s,f
   % lq(l,:,:)=0;
    for f=1:F
        expected_audio=lambda*circshift(A(:,f),round(alpha*l+beta));
        delta_audio=audio(:,f,2)-expected_audio;
        lq(l,s,f)=lq(l,s,f) -.5*asigma*(delta_audio'*delta_audio) + .5*sam*log(asigma/(2*pi));
        if ~O(f),
        y=frames(:,:,f);
            for s=1:S,
                mu_ls = mu_b.*(~circshift(gamma(:,:,s),[0,l]))...
                               + circshift(mu(:,:,s).*gamma(:,:,s),[0,l]);
                centered=y(:)-mu_ls(:);
   %         q(l,s,f)=q(l,s,f)*exp(-.5*vsigma*(centered'*centered));
                lq(l,s,f)=lq(l,s,f) + log(Pi(s)) - .5*vsigma*(centered'*centered) + .5*W*H*log(vsigma/(2*pi));
            end
        end
    end
end
if sum(isnan(q(:))) ~= 0,
    error('q contains NaN');
end


% Nomalization
q=reshape(lq,[W*S,F]);
a=max(q,[],1);
q=exp(bsxfun(@minus,q,a));
q=bsxfun(@rdivide,q,sum(q,1));
q=reshape(q,[W,S,F]);

% Nomalization
%q_partial=reshape(lq,[W,S,F]);
q_p = lq; % q with partial normalization
a=max(q_p,[],1);
q_p=exp(bsxfun(@minus,q_p,a));
q_p=bsxfun(@rdivide,q_p,sum(q_p,1)); %% CAREFUL not yet normalized in (l,s)
%q=reshape(q,[W,S,F]);


% imagesc(squeeze(q(:,1,:)+q(:,2,:)))
%% M-Step

%% mu_b
fprintf('Updating mu_b...')

mu_b=zeros(H,W); %numérateur de la formule (avant moyenne)
den=zeros(H,W); %dénominateur
for f=1:F
    if O(f), % occlusion
        mu_b=mu_b + frames(:,:,f);
        den=den + 1;
    else
        for s=1:S,
            i0=ftm2(q(:,s,f)',~gamma(:,:,s),1);
            mu_b=mu_b + i0.*frames(:,:,f);
            den=den + i0;
        end
    end
end
mu_b=mu_b ./ den;
fprintf(' Ok!\n')

%% Pi
fprintf('Updating Pi...')
Pi=mean(sum(q(:,:,~O),1),3);
if sum(Pi==0)~=0,
    Pi = ones(1,S)/S;
end
if all(isnan(Pi))
    error('Pi is all NaN...')
end
fprintf(' Ok!\n')

%% vsigma (Psi)
fprintf('Updating vsigma... :')
den=0;
for f=1:F
    y=frames(:,:,f);
    if O(f),
        cen=(y-mu_b);
        den=den+(cen(:)'*cen(:));
    else
        temp=y.^2; % développe le produit de e1 de manière vectorielle
        for s=1:S,
            temp=temp + mu_b.^2.*ftm2(q(:,s,f)',~gamma(:,:,s),1) ... %carres: ~gamma(:,:,s).^2=~gamma(:,:,s)
                + ftm2(q(:,s,f)',gamma(:,:,s).*mu(:,:,s).^2,1) ...
                - 2*y.*mu_b.*ftm2(q(:,s,f)',~gamma(:,:,s),1) ...%doubles produits, similarites avec ce qui precede
                - 2*y.*ftm2(q(:,s,f)',gamma(:,:,s).*mu(:,:,s),1); %
            % - 2*mu_b.*ftm2(q(:,s,f)'*ftm2(q(:,s,f)',~gamma(:,:,s),1);
            %              nul
        end
        %on transforme le vecteur en scalaire
        den=den+sum(temp(:));
    end
end
vsigma=H*W*F/den;
fprintf(' Ok!\n')

%% mu
fprintf('Updating mu... (depends on Pi)')
mu=zeros(H,W,S);

for f=1:F
    if ~O(f),
        for s=1:S,
            mu(:,:,s) = mu(:,:,s)+ftm2(q(:,s,f)',gamma(:,:,s),-1);
        end
    end
end
for s=1:S,
    mu(:,:,s)=mu(:,:,s)/(sum(~O)*Pi(s));
end
% already normalized accordingly using q_p instead of q.
fprintf(' Ok!\n')


%% O
fprintf('Updating O...')
e0=zeros(F,1);
e1=zeros(F,1);
for f=1:F
    y=frames(:,:,f);
    cen=(y-mu_b);   
    e0(f)=(cen(:)'*cen(:));
    temp=y.^2; % développe le produit de e1 de manière vectorielle
    for s=1:S,
        temp=temp + mu_b.^2.*ftm2(q(:,s,f)',~gamma(:,:,s),1) ... %carres: ~gamma(:,:,s).^2=~gamma(:,:,s)
            + ftm2(q(:,s,f)',gamma(:,:,s).*mu(:,:,s).^2,1) ...
            - 2*y.*mu_b.*ftm2(q(:,s,f)',~gamma(:,:,s),1) ...%doubles produits, similarites avec ce qui precede
            - 2*y.*ftm2(q(:,s,f)',gamma(:,:,s).*mu(:,:,s),1); %
        % - 2*mu_b.*ftm2(q(:,s,f)'*ftm2(q(:,s,f)',~gamma(:,:,s),1);
        %              nul
    end
    %on transforme le vecteur en scalaire
    e1(f)=sum(temp(:));
end
O=e0>e1;
fprintf(' Ok!\n')

%% Gamma
gamma_hist(:,:,:,MASTER)=gamma;

fprintf('Updating gamma...')
for s=1:S
    e0=zeros(H,W);
    e1=zeros(H,W);
    for f=1:F,
    y=frames(:,:,f);
        if ~O(f),
            e0=e0+ftm2(q(:,s,f)',(frames(:,:,f)-mu_b).^2,-1);
            e1=e1 + ftm2(q(:,s,f)',(frames(:,:,f)).^2,-1) ...
                  + sum(q(:,s,f))*mu(:,:,s).^2 ...
                 -2*ftm2(q(:,s,f)',(frames(:,:,f)),-1).*mu(:,:,s);
        end
    end
%gamma(:,:,s)=e0>e1
gamma(:,:,s)=e0>e1;
end
fprintf(' Ok!\n')

if all(gamma == 0),
    fprintf('gamma is 0. Resetting to initial gamma...\n'); 
    gamma=gamma_hist(:,:,:,1);
end

%% %%%%%%%%%%% AUDIO %%%%%%%%%%%%%%%%

%% Recentrer gammas, mu and q (nécessaire pour avoir une correspondance en l)
[ gamma, mu, q ] = recenter( gamma, mu, q, 10 );

q_l=sum(q,2); % ne dépend plus de s
% dilataion de q(l) pour avoir q(s)
% où s est le décalage en samples (changement d'unités)
q_s=zeros(1,sam,F);
for f=1:F,
    q_s(1,:,f)=interp1(1:W,q_l(:,1,f),((1:sam)'-beta)/alpha,'linear',0);
    q_s(1,:,f)=q_s(1,:,f)/sum(q_s(1,:,f)); % renormalization dans la nouvelle échelle
end

%% Eta
fprintf('Updating Eta...')
den=0;
for f=1:F
    yaudio = lambda*A(:,f);
        temp = yaudio.^2; % développe le produit de e1 de manière vectorielle
        for s=1:S,
            temp = temp + ftm2(q_s(1,:,f),audio(:,f,2)'.^2,-1)' ...
                - 2*yaudio.*ftm2(q_s(1,:,f),audio(:,f,2)',-1)';
        end
        %on transforme le vecteur en scalaire
        den=den+sum(temp(:));
end
%asigma=sam*F/den;
fprintf(' Disabled! asigma=%f\n',asigma)

%% lambda
fprintf('Updating lambda...')
%lambda=0;
%for f=1:F
%    lambda=lambda+( ftm2(q_s(1,:,f),audio(:,f,2)',-1) * A(:,f)...
%                    / (A(:,f)'*A(:,f)) );
%end
%lambda=lambda/F;
fprintf(' Disabled! Ok!\n')

%% A_f
fprintf('Updating A_f...')
%for f=1:F
%    A(:,f)=(audio(:,f,1) + lambda * ftm2(q_s(1,:,f), audio(:,f,2)', -1)')...
%           / (1+lambda^2);
%end
fprintf(' Disabled! Ok!\n')

%% alpha and beta
fprintf('Updating alpha and beta...')
best_l=zeros(F,1);
% position la plus probable trouvée pour chaque frame
for f=1:F,
    [~,best_l(f)]=max(sum(q(:,:,f),2),[],1);
end
cleaning = abs(sau)<100;
clean_best_l = best_l(cleaning);
clean_sau = sau(cleaning);
[ p, inliers ] = ransacline( [clean_best_l'; clean_sau'], 100, 10 );
%alpha=p(1);
%beta=p(2);
fprintf('Disabled! Ok!\n')
fprintf('alpha: %f, beta: %f\n', alpha, beta)

imshow(Gamma)
drawnow
end
figure
imagesc(Qs)
%imshow(cat(3,squeeze(q(:,1,:)),squeeze(q(:,2,:)),zeros(W,F)))