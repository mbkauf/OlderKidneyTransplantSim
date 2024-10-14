get_psa_qalys <- function(num) {
    ## Hanmer 2006 parameters
    # Male (alphas and betas)
    m_40_alpha <- 0.887*(((0.887*(1-0.887))/0.000013)-1)
    m_40_beta  <- (1-0.887) * (((0.887*(1-0.887))/0.000013)-1)
    m_60_alpha <- 0.840*(((0.840*(1-0.840))/0.000037)-1)
    m_60_beta  <- (1-0.840) * (((0.840*(1-0.840))/0.000037)-1)
    m_70_alpha <- 0.802*(((0.802*(1-0.802))/0.000051)-1)
    m_70_beta  <- (1-0.802) * (((0.802*(1-0.802))/0.000051)-1)
    m_80_alpha <- 0.782*(((0.782*(1-0.782))/0.000163)-1)
    m_80_beta  <- (1-0.782) * (((0.782*(1-0.782))/0.000163)-1)

    # Female (alphas and betas)
    f_40_alpha <- 0.863*(((0.863*(1-0.863))/0.000017)-1)
    f_40_beta  <- (1-0.863) * (((0.863*(1-0.863))/0.000017)-1)
    f_60_alpha <- 0.811*(((0.811*(1-0.811))/0.000031)-1)
    f_60_beta  <- (1-0.811) * (((0.811*(1-0.811))/0.000031)-1)
    f_70_alpha <- 0.771*(((0.771*(1-0.771))/0.000044)-1)
    f_70_beta  <- (1-0.771) * (((0.771*(1-0.771))/0.000044)-1)
    f_80_alpha <- 0.724*(((0.724*(1-0.724))/0.000138)-1)
    f_80_beta  <- (1-0.724) * (((0.724*(1-0.724))/0.000138)-1)

    # Get vector of QALY values for Hanmer paper
    samp_m    <- ordered_qalys(n.runs = num,
                               alphas = c(m_40_alpha, m_60_alpha,
                                          m_70_alpha, m_80_alpha),
                               betas = c(m_40_beta, m_60_beta,
                                         m_70_beta, m_80_beta))
    samp_f    <- ordered_qalys(n.runs = num,
                               alphas = c(f_40_alpha, f_60_alpha,
                                          f_70_alpha, f_80_alpha),
                               betas = c(f_40_beta, f_60_beta,
                                         f_70_beta, f_80_beta))
    samp_m_40 <- samp_m[, 1]
    samp_m_60 <- samp_m[, 2]
    samp_m_70 <- samp_m[, 3]
    samp_m_80 <- samp_m[, 4]
    samp_f_40 <- samp_f[, 1]
    samp_f_60 <- samp_f[, 2]
    samp_f_70 <- samp_f[, 3]
    samp_f_80 <- samp_f[, 4]

    ## Wyld 2012 parameters
    # Sample sizes
    n_pretx <- 343
    n_03    <- 598
    n_48    <- 190
    n_912   <- 686
    n_13    <- 509

    # alphas and betas
    pretx_alpha    <- 0.50 * n_pretx
    pretx_beta     <- (1 - 0.50) * n_pretx
    post_03_alpha  <- 0.63 * n_03
    post_03_beta   <- (1 - 0.63) * n_03
    post_48_alpha  <- 0.71 * n_48
    post_48_beta   <- (1 - 0.71) * n_48
    post_912_alpha <- 0.64 * n_912
    post_912_beta  <- (1 - 0.64) * n_912
    post_13_alpha  <- 0.63 * n_13
    post_13_beta   <- (1 - 0.63) * n_13

    # Get vector of QALY values for Wyld papers
    samp_pretx    <- rbeta(num, pretx_alpha, pretx_beta)
    samp_post_03  <- rbeta(num, post_03_alpha, post_03_beta)
    samp_post_48  <- rbeta(num, post_48_alpha, post_48_beta)
    samp_post_912 <- rbeta(num, post_912_alpha, post_912_beta)
    samp_post_13  <- rbeta(num, post_13_alpha, post_13_beta)

    ## Begin calculations for QALYs
    diff_m_pretx    <- samp_pretx - samp_m_40
    diff_m_post_03  <- samp_post_03 - samp_m_40
    diff_m_post_48  <- samp_post_48 - samp_m_40
    diff_m_post_912 <- samp_post_912 - samp_m_40
    diff_m_post_13  <- samp_post_13 - samp_m_40

    diff_f_pretx    <- samp_pretx - samp_f_40
    diff_f_post_03  <- samp_post_03 - samp_f_40
    diff_f_post_48  <- samp_post_48 - samp_f_40
    diff_f_post_912 <- samp_post_912 - samp_f_40
    diff_f_post_13  <- samp_post_13 - samp_f_40

    u_dialysis_m_60   <- samp_m_60 + diff_m_pretx
    u_posttx_03_m_60  <- samp_m_60 + diff_m_post_03
    u_posttx_48_m_60  <- samp_m_60 + diff_m_post_48
    u_posttx_912_m_60 <- samp_m_60 + diff_m_post_912
    u_posttx_13_m_60  <- samp_m_60 + diff_m_post_13

    u_dialysis_f_60   <- samp_f_60 + diff_f_pretx
    u_posttx_03_f_60  <- samp_f_60 + diff_f_post_03
    u_posttx_48_f_60  <- samp_f_60 + diff_f_post_48
    u_posttx_912_f_60 <- samp_f_60 + diff_f_post_912
    u_posttx_13_f_60  <- samp_f_60 + diff_f_post_13

    u_dialysis_m_70   <- samp_m_70 + diff_m_pretx
    u_posttx_03_m_70  <- samp_m_70 + diff_m_post_03
    u_posttx_48_m_70  <- samp_m_70 + diff_m_post_48
    u_posttx_912_m_70 <- samp_m_70 + diff_m_post_912
    u_posttx_13_m_70  <- samp_m_70 + diff_m_post_13

    u_dialysis_f_70   <- samp_f_70 + diff_f_pretx
    u_posttx_03_f_70  <- samp_f_70 + diff_f_post_03
    u_posttx_48_f_70  <- samp_f_70 + diff_f_post_48
    u_posttx_912_f_70 <- samp_f_70 + diff_f_post_912
    u_posttx_13_f_70  <- samp_f_70 + diff_f_post_13

    u_dialysis_m_80   <- samp_m_80 + diff_m_pretx
    u_posttx_03_m_80  <- samp_m_80 + diff_m_post_03
    u_posttx_48_m_80  <- samp_m_80 + diff_m_post_48
    u_posttx_912_m_80 <- samp_m_80 + diff_m_post_912
    u_posttx_13_m_80  <- samp_m_80 + diff_m_post_13

    u_dialysis_f_80   <- samp_f_80 + diff_f_pretx
    u_posttx_03_f_80  <- samp_f_80 + diff_f_post_03
    u_posttx_48_f_80  <- samp_f_80 + diff_f_post_48
    u_posttx_912_f_80 <- samp_f_80 + diff_f_post_912
    u_posttx_13_f_80  <- samp_f_80 + diff_f_post_13

    v_u_dialysis_m <- cbind(u_dialysis_m_60, u_dialysis_m_70, u_dialysis_m_80)
    v_u_dialysis_f <- cbind(u_dialysis_f_60, u_dialysis_f_70, u_dialysis_f_80)
    v_u_posttx_03_m <- cbind(u_posttx_03_m_60, u_posttx_03_m_70, u_posttx_03_m_80)
    v_u_posttx_03_f <- cbind(u_posttx_03_f_60, u_posttx_03_f_70, u_posttx_03_f_80)
    v_u_posttx_48_m <- cbind(u_posttx_48_m_60, u_posttx_48_m_70, u_posttx_48_m_80)
    v_u_posttx_48_f <- cbind(u_posttx_48_f_60, u_posttx_48_f_70, u_posttx_48_f_80)
    v_u_posttx_912_m <- cbind(u_posttx_912_m_60, u_posttx_912_m_70, u_posttx_912_m_80)
    v_u_posttx_912_f <- cbind(u_posttx_912_f_60, u_posttx_912_f_70, u_posttx_912_f_80)
    v_u_posttx_13_m <- cbind(u_posttx_13_m_60, u_posttx_13_m_70, u_posttx_13_m_80)
    v_u_posttx_13_f <- cbind(u_posttx_13_f_60, u_posttx_13_f_70, u_posttx_13_f_80)

    return(list(v_u_dialysis_m, v_u_dialysis_f, v_u_posttx_03_m, v_u_posttx_03_f,
           v_u_posttx_48_m, v_u_posttx_48_f, v_u_posttx_912_m, v_u_posttx_912_f,
           v_u_posttx_13_m, v_u_posttx_13_f))
}

l_qalys <- get_psa_qalys(10000)
all(l_qalys[[1]] < l_qalys[[3]])
all(l_qalys[[1]] < l_qalys[[5]])
all(l_qalys[[1]] < l_qalys[[7]])
all(l_qalys[[1]] < l_qalys[[9]])

all(l_qalys[[2]] < l_qalys[[4]])
all(l_qalys[[2]] < l_qalys[[6]])
all(l_qalys[[2]] < l_qalys[[8]])
all(l_qalys[[2]] < l_qalys[[10]])
