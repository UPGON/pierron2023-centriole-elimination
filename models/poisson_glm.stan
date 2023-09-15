data {
        int<lower=1> N;
        int<lower=1> J;
        vector[N] x;
        array[N] int<lower=0>  y;
}

parameters {
    real b0;
    real b1;
}

model {
    b0 ~ normal(0, 10);
    b1 ~ normal(0, 10);
    y ~ poisson_log(b0 + b1 * x);
}
