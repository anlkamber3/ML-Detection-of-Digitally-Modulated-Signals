clear all
warning off
%%% BPSK over AWGN
%%%%% SIGNAL CONSTELLATION %%%%%%%%%
symbolBook=[1 1i];
bitBook=[0;1];
nBitPerSym=size(bitBook,2);
M=length(symbolBook);
%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
nSymPerFrame=3000;
nBitsPerFrame=nSymPerFrame*nBitPerSym;
max_nFrame=2000;
fErrLim=100;
snr_db=0:20;
nBitErrors=zeros(length(snr_db), 1);
nTransmittedFrames=zeros(length(snr_db), 1);
nErroneusFrames=zeros(length(snr_db), 1);
SYMBOLBOOK=repmat(transpose(symbolBook),1,nSymPerFrame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nEN = 1:length(snr_db) % SNR POINTS
    this_snr=snr_db(nEN);
    sigma_noise = 1/sqrt(10^(this_snr/10));
    while (nTransmittedFrames(nEN)<max_nFrame) && (nErroneusFrames(nEN)<fErrLim)
        nTransmittedFrames(nEN) = nTransmittedFrames(nEN) + 1;
        %%%%%%%%%% INFORMATION GENERATION %%%%%%%%%%
        trSymIndices=randi(M,[1,nSymPerFrame]);
        trSymVec=symbolBook(trSymIndices);
        trBitsMat=bitBook(trSymIndices,:)';
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        noise=1/sqrt(2)*[randn(1, nSymPerFrame) + 1j*randn(1,nSymPerFrame)];
        recSigVec=trSymVec+sigma_noise*noise;
        %%%% DETECTOR %%%%%%%%%%%%
        RECSIGVEC=repmat(recSigVec,length(symbolBook),1);
        distance_mat=abs(SYMBOLBOOK-RECSIGVEC);
        [~, det_sym_ind]=min(distance_mat,[],1);
        detected_bits=[bitBook(det_sym_ind, :)]';
        err = sum(sum(abs(trBitsMat-detected_bits)));
        nBitErrors(nEN)=nBitErrors(nEN)+err;
        if err~=0
            nErroneusFrames(nEN)=nErroneusFrames(nEN)+1;
        end
    end % End of while loop
end %end for (SNR points)
semilogy(snr_db, nBitErrors./nTransmittedFrames/nBitsPerFrame,"r-x");
legend("BPSK BER Performance versus SNR",'interpreter','latex');
grid("on");
xlabel('SNR (dB)'); ylabel('BER');
