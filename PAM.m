clear all
warning off
%%% 4-PAM over AWGN
%%%%% SIGNAL CONSTELLATION %%%%%%%%%
symbolBook=[-3 -1 1 3]; %To make average energy unity, multiplier equals to sqrt(1/5)
symbolBook = sqrt(1/5) * symbolBook;

bitBook_uniform= [0 0;0 1;1 0;1 1];
bitBook_grey = [0 0;0 1;1 1;1 0];

nBitPerSym=size(bitBook_uniform,2);

M=length(symbolBook);

%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
nSymPerFrame=3000;
nBitsPerFrame=nSymPerFrame*nBitPerSym;
max_nFrame=2000;
fErrLim=100;
snr_db=0:20;

%%%%%%%%%%%%%%%%%%%%%%GREY MAPPING%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBitErrors_grey=zeros(length(snr_db), 1);
nTransmittedFrames_grey=zeros(length(snr_db), 1);
nErroneusFrames_grey=zeros(length(snr_db), 1);

%%%%%%%%%%%%%%%%%%%%%%UNIFORM MAPPING%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBitErrors_uniform=zeros(length(snr_db), 1);
nTransmittedFrames_uniform=zeros(length(snr_db), 1);
nErroneusFrames_uniform=zeros(length(snr_db), 1);

SYMBOLBOOK=repmat(transpose(symbolBook),1,nSymPerFrame);


for nEN = 1:length(snr_db) % SNR POINTS
    this_snr=snr_db(nEN);
    sigma_noise = 1/sqrt(10^(this_snr/10));
    %%%%%%%%%%%%%%%%%%%%GREY MAPPING%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (nTransmittedFrames_grey(nEN)<max_nFrame) && (nErroneusFrames_grey(nEN)<fErrLim)
        nTransmittedFrames_grey(nEN) = nTransmittedFrames_grey(nEN) + 1;
        %%%%%%%%%% INFORMATION GENERATION %%%%%%%%%%
        trSymIndices=randi(M,[1,nSymPerFrame]);
        trSymVec=symbolBook(trSymIndices);
        trBitsMat=bitBook_grey(trSymIndices,:)';
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        noise=1/sqrt(2)*[randn(1, nSymPerFrame) + 1j*randn(1,nSymPerFrame)];
        recSigVec=trSymVec+sigma_noise*noise;
        %%%% DETECTOR %%%%%%%%%%%%
        RECSIGVEC=repmat(recSigVec,length(symbolBook),1);
        distance_mat=abs(SYMBOLBOOK-RECSIGVEC);
        [~, det_sym_ind]=min(distance_mat,[],1);
        detected_bits=[bitBook_grey(det_sym_ind, :)]';
        err= sum(sum(abs(trBitsMat-detected_bits)));
        nBitErrors_grey(nEN)=nBitErrors_grey(nEN)+err;

        if err~=0
            nErroneusFrames_grey(nEN)=nErroneusFrames_grey(nEN)+1;
        end
    end % End of while loop
    %%%%%%%%%%%%%%%%%%%%%%%%%UNIFORM MAPPING%%%%%%%%%%%%%%%%%%%%%%%
    while (nTransmittedFrames_uniform(nEN)<max_nFrame) && (nErroneusFrames_uniform(nEN)<fErrLim)
        nTransmittedFrames_uniform(nEN) = nTransmittedFrames_uniform(nEN) + 1;
        %%%%%%%%%% INFORMATION GENERATION %%%%%%%%%%
        trSymIndices=randi(M,[1,nSymPerFrame]);
        trSymVec=symbolBook(trSymIndices);
        trBitsMat=bitBook_uniform(trSymIndices,:)';
        %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
        noise=1/sqrt(2)*[randn(1, nSymPerFrame) + 1j*randn(1,nSymPerFrame)];
        recSigVec=trSymVec+sigma_noise*noise;
        %%%% DETECTOR %%%%%%%%%%%%
        RECSIGVEC=repmat(recSigVec,length(symbolBook),1);
        distance_mat=abs(SYMBOLBOOK-RECSIGVEC);
        [~, det_sym_ind]=min(distance_mat,[],1);
        detected_bits=[bitBook_uniform(det_sym_ind, :)]';
        err = sum(sum(abs(trBitsMat-detected_bits)));
        nBitErrors_uniform(nEN)=nBitErrors_uniform(nEN)+err;
        if err~=0 
            nErroneusFrames_uniform(nEN)=nErroneusFrames_uniform(nEN)+1;
        end
    end % End of while loop
end %end for (SNR points)

snr = 10.^(snr_db/10);
sim_res=[nBitErrors_uniform nTransmittedFrames_uniform;nBitErrors_grey nTransmittedFrames_grey];
theory_uniform = qfunc(sqrt(6/15).*sqrt(snr));
theory_grey = (3/4)*qfunc(sqrt((6/15)*snr));
semilogy(snr_db, nBitErrors_grey./nTransmittedFrames_grey/nBitsPerFrame, 'b-x',snr_db, nBitErrors_uniform./nTransmittedFrames_uniform/nBitsPerFrame,"r-x",snr_db,theory_uniform,"rs",snr_db,theory_grey,"bd");
legend("4-PAM Grey Mapping","4-PAM Uniform Mapping","Theory of Uniform Mapping","Theory of Grey Mapping",'interpreter','latex');
grid("on");
xlabel('SNR (dB)'); ylabel('BER');
