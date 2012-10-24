function Z = Bh_to_Z(Bh)

    'impossible'
    die

  N = length(Bh);
  Extr = zeros(N,2*N);
  for j = 1:N
    Extr(j,j) = 1;
  end

  Z1 = fft(Bh);
  M = ComputeM(N);
  Lambda = ComputeLambda(nbpoints_discr,step,H);
  Z2 = Lambda^(-1)*Z1;
  Z3 = M^(-1) * Z2;
  B1 = fft(Bh);
  Lambda = ComputeLambda(nbpoints_discr,step,H);
  Z = (Lambda*M)^(-1)*B1; %code efficiently

  