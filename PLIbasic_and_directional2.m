function [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(ffteeg, N, F, overlap)
% Calculates coherence, coherency, and phase locking measures
% 
% Input :   ffteeg - complex FFT, Nsegments x Nchannels x Nfrequencies
%           N - number of segments to use, picked up at random. Should be less or equal to size(ffteeg,1)
%           F - frequencies to use, e.g.  [1:50], in points, not Hz !
%           overlap - Nsegments x Nsegments matrix with '0' indicating no
%                     overlap and '1'  indicating overlap between segments; if no
%                     overlap - eye(Nsegments); if [] -   no dbWPLI,  ubPLI calculated

% Output:   Coh - coherence (0 to 1)
%           Coherency -  complex coherency
%           CoherencyZ - normalized coherency (by Nolte et al, 2004)
%           PLI, WPLI, ubPLI, dbWPLI - phase locking measures
%
% Copyright (C) 2013 Elena Orekhova, Birkbeck, University of London

% Adapted by Cosmin Stamate and Rianne Haartsen, 2017, Birkbeck, University
% London
%
    [ll,kk] = ind2sub([N,N],1:N^2);
    isBelowDiagonal = ll > kk;
    overlap(isBelowDiagonal) = 1;
    
    B = repmat(overlap,[1,1, size(F,2)]);


    if isempty(overlap)==1
        segments = randperm(size(ffteeg,1));     
        segments = segments(1:N); 
        dbWPLI = [];
    else
        segments = 1:N; 
        dbWPLI = zeros(size(F,2), size(ffteeg,2),size(ffteeg,2));
    end
    
    ffteeg = ffteeg(segments,:,F);
    xy = zeros( N, size(ffteeg,3),size(ffteeg,2),size(ffteeg,2) ); % trials x Fs x channels x channels


   
    for i=1:size(ffteeg,2)   %channels
        i
        
       for ii=1:size(ffteeg,2)  %channels
          
           tmp1 = zeros(size (xy,1),size (xy,1), size(xy,2)); % N x N x freq 
           tmp2 = zeros(size (xy,1),size (xy,1), size(xy,2)); % N x N x freq 
           xy(:,:,i,ii)=squeeze(ffteeg(:,i,:)).*conj(squeeze(ffteeg(:,ii,:))); % N x Fs x ch x ch
           
           if isempty(overlap)==0 % Calculate dbWPLI and  ubPLI only if have info about segments overlap!
               
               for j=1:N  %epochs                   
                       a = repmat(imag(xy(j,:,i,ii)), N, 1);
                       b = imag(xy(:,:,i,ii));
                       tmp1(j,:,:)=a.*b;
                       tmp2(j,:,:)=sign(a).*sign(b);                      
               end
              
              tmp1(B==1) = 0;  tmp2(B==1) = 0; % exclude overlapping segments 
              
              tmp1_abs = abs(tmp1);

              dbWPLI(:,i,ii)= sum(squeeze(sum(tmp1,1)),1)./sum(squeeze(sum(tmp1_abs,1)),1); % debiased  WPLI estimator % modified from  Vinck et al. (2011)
              ubPLI(:,i,ii) = sum(squeeze(sum(tmp2,1)),1)./size(find(tmp2(:,:,1)),1);
           end
           
          PLI(:,i,ii) = shiftdim(  abs(mean(sign(imag(xy(:,:,i,ii))),1))    );  % Stam et al. (2007)
          WPLI(:,i,ii)= shiftdim(  abs(     mean(imag(xy(:,:,i,ii))))  ./  mean (abs(imag(xy(:,:,i,ii))))      );  % Vinck et al. (2011)
       end  
       
    end
    
    WPLI(find(isnan(WPLI)))=0; 
    dbWPLI(find(isnan(dbWPLI)))=0;
    PLI(find(isnan(PLI)))=0; 
    ubPLI(find(isnan(ubPLI)))=0;
    
    ubPLI(:, find(eye(size(ubPLI,3))) )=0; % modified unbiased squared PLI estimator % Vinck et al. (2011, 2012)

end

