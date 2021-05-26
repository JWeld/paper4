#STAN GLMS####
library(rstanarm)
library(bayesplot)
library(broom)
library(tidybayes)
library(modelr)
library(sjPlot)
library(parallel)

dat <- all_data2
#N


glm1 <- stan_glmer(pc_dist_base ~ NH4M * NO3M + SO4SM + PREC + latitude +
                       longitude + survey_year + (1 |ID_site), 
                       chains = 4, cores = 4, iter= 6000, data = dat, adapt_delta = 0.95)

glm2 <- stan_glmer(pc_dist_base ~ NTOT + SO4SM + PREC + latitude +
                     longitude + survey_year + (1 |ID_site), 
                   chains = 4, cores = 4, iter= 6000, data = dat, adapt_delta = 0.95)

glm2.5 <- stan_glmer(dispersion ~ NH4M * NO3M + SO4SM + TEMP + latitude + longitude +
                     PREC + survey_year + (1|ID_site), 
                   chains = 4, cores = 4, iter= 6000, data = dat, adapt_delta = 0.95)

#summaries stan
fit <- glm2
launch_shinystan(fit)
#model comparison
fit1 <- glm1
fit2 <- glm2

loo1 <- loo(fit1, cores = 6)
print(loo1)
loo2 <- loo(fit2, cores = 6)
print(loo2)


plot(loo1, label_points = TRUE)
plot(loo2, label_points = TRUE)
loo_compare(loo1, loo2)

#fit <- stan_glms_job_results$stan_glmer2.5.br
summary(fit2)
prior_summary(fit2)
#library(MuMIn)
r.squaredGLMM(fit2) #MuMin

#plot fitted vs data on map
plot(x,y,cex=0.1*fitted(fit),col="gray",asp=1)
plot(x,y,cex=0.1*dat$rich.y,col="blue",asp=1)

rsq <- bayes_R2(fit)
print(median(rsq))

summary(fit, 
        pars = c("(Intercept)", "sigma", "Sigma[grp_tree_species:(Intercept),(Intercept)]"),
        probs = c(0.025, 0.975),
        digits = 2)

nd <- select(simple.FD,R.y, NH4M, NO3M, plot.x, ID_siteonly, country, survey_year)
ytilde <- posterior_predict(fit, nd, draws = 500)

print(dim(ytilde))  # 500 by 3 matrix (draws by nrow(nd))
ytilde <- data.frame(count = c(ytilde),
                     outcome = rep(simple.FD$R.y, each = 500))

ggplot2::ggplot(ytilde, ggplot2::aes(x=outcome, y=count)) +
  ggplot2::geom_point() +
  ggplot2::ylab("predicted count")

plot(ytilde$count,ytilde$outcome)
       
#To use the posterior draws with the functions in the bayesplot package
#we’ll extract them from the fitted model object.We’ve used as.array above
#(as opposed to as.matrix) because it keeps the Markov chains separate 
#(stan_glm runs four chains by default). 

posterior2 <- as.array(fit)


mcmc_intervals(posterior2, pars = c("NH4M", "(Intercept)"))

mcmc_intervals(posterior2, pars = c("NO3M" , "NH4M", "(Intercept)", "NH4M:NO3M", "survey_year" ))
mcmc_hex(posterior2, pars = c("NO3M" , "NH4M"))
mcmc_areas(posterior2, 
           pars = c("NO3M" , "NH4M", "NH4M:NO3M"),
           prob = 0.8, # 80% intervals
           prob_outer = 0.99, # 99%
           point_est = "mean"
)

plot(fit)
plot(fit, "trace")
pp_check(fit,plotfun = "stat", stat = "mean")
pp_check(fit, nreps = 100)

#plot model predictions and actuals data

(model_fit <- dat %>%
    data_grid(NH4M = seq_range(NH4M, n = 101), NO3M = seq_range(NO3M, n = 101)) %>%
    add_predicted_draws(fit) %>%
    ggplot(aes(x = NH4M, y = dispersion)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),
                    alpha = 1/2, colour = "black")
  +
    geom_point(data = simple.FD, colour = "darkseagreen4", size = 3) +
    scale_fill_brewer(palette = "Greys"))


launch_shinystan(fit)

y_pred = posterior_predict(fit)

#summary from stan object
bh_summary <- fit$stan_summary %>% 
  as_tibble() %>% 
  mutate(variable = rownames(.)) %>% 
  select(variable, everything()) %>% 
  as_tibble()

bh_summary %>% head()


bh_summary %>% 
  ggplot(aes(n_eff)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = 4000), color = 'red') #4 chains, with 2000 iterations, 
#half of which are warmup, meaning we sample 1000 iterations in each chain, 
#so the max n_eff possible in this case is 4000

bh_summary %>% 
  ggplot(aes(Rhat)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = 1.1), color = 'red') #Gelman recommends that Rhat for each parameter be less than 1.1

bh_summary %>% 
  filter(variable %in% c('(Intercept)','NH4M','NO3M')) %>% 
  ggplot() + 
  geom_linerange(aes(variable, ymin = `2.5%`,ymax = `97.5%`)) + 
  geom_crossbar(aes(variable, mean, ymin = `25%`, ymax = `75%`), fill= 'grey') + 
  facet_wrap(~variable, scales = 'free')



# pvp <- posterior_vs_prior(fit, color_by = "vs", group_by = TRUE)
# 
#   pvp + 
#   ggplot2::geom_hline(yintercept = 0, size = 0.3, linetype = 3) + 
#   ggplot2::coord_flip() + 
#   ggplot2::ggtitle("Comparing the prior and posterior")

# show group-level (varying) parameters and group by parameter
posterior_vs_prior(pvp, pars = "varying",
                   group_by_parameter = TRUE, color_by = "vs")

# group by parameter and allow axis scales to vary across facets
posterior_vs_prior(pvp, regex_pars = "period",
                   group_by_parameter = TRUE, color_by = "none",
                   facet_args = list(scales = "free"))


newdata = dat %>% cbind(t(y_pred)) %>% gather(key = "Rep", value = "Value",
                                               -NH4M, -rich.y)

ggplot(dat, aes(Value, x = N.y)) + geom_violin(color = "blue",
                                                fill = "blue", alpha = 0.5) + geom_violin(data = dat1, aes(y = NH4M,
                                                                                                           x = N.y), fill = "red", color = "red", alpha = 0.5)


######plotting regression lines ##
# Extract the (post-warmup) posterior draws
posterior1 <- as.matrix(fit)
colnames(posterior1)
means1 <- colMeans(posterior1)

# Take random 100 posterior draws of intercept and slope
# 100 isn't special, but enough to show uncertainty without
# making the graph unreadable
betas <- posterior1[sample(nrow(posterior1), 200), 1:2]


# Plot regression lines implied by the betas
blues <- color_scheme_get("brightblue")
mod1p1 <- ggplot(dat, aes(x = NH4M, y = dispersion)) +
  geom_point(color = "gray30") +
  geom_abline(
    intercept = betas[, 1], 
    slope = betas[, 2],
    color = blues[[2]], 
    size = 0.15, 
    alpha = 0.5
  ) +
  geom_abline(
    intercept = means1[1], 
    slope = means1[2],
    size = 1.25, 
    color = blues[[6]]
  ) +
  ylim(0, 10)

plot(mod1p1)
#ggsave("plots/regression1.pdf", width = 8, height = 6)
temp <- drop_na(dat)
(model_fit <- simple.FD %>%
    data_grid(nh4 = NH4M) %>%
    add_predicted_draws(fit) %>%
    ggplot(aes(x = NH4M, y = rich.y)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),
                    alpha = 1/2, colour = "black") +
    geom_point(data = temp, colour = "darkseagreen4", size = 3) +
    scale_fill_brewer(palette = "Greys"))


# Optional 
#    Estimate out-of-sample predictive performance
#    Can use later to compare to the hierarchical model
looSingle <- loo(SingleLevelModel) 