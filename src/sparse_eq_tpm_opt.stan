data {
  int TR;// number of TRs
  int R;// number of reads
  int V;
  int Ulen;
  vector [V] w ;
  int  v [V];
  int   u [Ulen];
  vector[TR] trlen;
  vector [R] classweights;
}

parameters {
  simplex [TR] tpm; // the ratio of steady state ratio to ribo
}

transformed parameters {
  vector [TR] n;
  n = tpm ./ trlen;
}

model {
    target += log(
      csr_matrix_times_vector(R, TR, w, v, u, n)
      // (map * n)
    ).* classweights;
}
