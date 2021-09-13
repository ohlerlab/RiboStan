data {
  int TR;// number of TRs
  int R;// number of reads
  matrix[R,TR] map;
  vector[TR] trlen;
  matrix [R,1] classweights;
}

parameters {
  simplex [TR] tpm; // the ratio of steady state ratio to ribo
}

transformed parameters {
  matrix [TR,1] n;
  n[,1] = tpm ./ trlen;
}

model {
    target += log(
      (map * n)
      // (map * n)
    ).* classweights;
}
