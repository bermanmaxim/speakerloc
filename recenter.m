function [ gamma, mu, q ] = recenter( gamma, mu, q, w )
%RECENTER Aligner les masques au d�but de l'image pour calculer l.
% w pr�cise la taille de la fen�tre utilis� pour le filtre de d�tection
% de saut.
for s=1:size(gamma,3),
    hpos = 2*mean(gamma(:,:,s),1)-1;
    %hpos est une version 1D de gamma qui varie entre -1 et 1
    filtre = [ones(1,w) -ones(1,w)]; % on veut d�tecter un �chelon positif
    result = cconv(hpos,filtre,size(hpos,2));
    [~, orig] = max(result);
    orig = orig - w;
    gamma(:,:,s)=circshift(gamma(:,:,s),[0 -orig]);
    mu(:,:,s)=circshift(mu(:,:,s),[0 -orig]);
    q(:,s,:)=circshift(q(:,s,:),[-orig 0 0]);
end

end

