function model = motif2Model(motif)
    A = [1 0 0 0];
    C = [0 1 0 0];
    G = [0 0 1 0];
    T = [0 0 0 1];
    M = [1 1 0 0];
    R = [1 0 1 0];
    W = [1 0 0 1];
    S = [0 1 1 0];
    Y = [0 1 0 1];
    K = [0 0 1 1];
    V = [1 1 1 0];
    H = [1 1 0 1];
    D = [1 0 1 1];
    B = [0 1 1 1];
    N = [1 1 1 1];
    X = [1 1 1 1];
    
    alphabet = ['ACGTMRWSYKVHDBNX'];
    matrices = [A;C;G;T;M;R;W;S;Y;K;V;H;D;B;N;X];
    for i=1:numel(motif)
        j = find(alphabet == upper(motif(i)));
        model.emat(i,:) = 1-matrices(j,:);
    end
    model.cutoff = 1;
end