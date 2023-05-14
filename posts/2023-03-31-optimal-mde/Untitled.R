library(tidyverse)
library(pwr)


samp_size <- Vectorize(function(mde, baseline, n_uuids, max_exp_per_year){
  
  # Group size which ensures 24 experiments per year
  n_min = n_uuids/max_exp_per_year/2
  # Estimated sample size for the mde
  n_est = ceiling(pwr.2p.test(h = ES.h(mde*baseline, baseline), power = 0.8)$n)
  # Take the max.  If n_est is smaller than n_min, then we would just run the experiment longer since
  # we can't run more of them.
  pmax( n_min, n_est )
})

one_sided_power = function(real_lift, n_per_group){
  # Only interested in the case when the estimated lift is
  # Greater than one, which corresponds to a one sided test.
  # However, you always run 2 tailed tests, so the significance level
  # is half of what is typically is.
  pwr.2p.test(h = ES.h(real_lift*baseline, baseline), 
              n = n_per_group,
              alternative = 'greater',
              sig.level = 0.025
  )$power
}


simr <- Vectorize(function(n_experiments_per_year, n_per_group, latent_lift, latent_sd){
   
  simulations <- crossing(
    sim = 1:1000, 
    experiment = 1:n_experiments_per_year,
    n_per_group = n_per_group
  )
  
  simulations %>% 
    mutate(
      # draw a real lift for each experiment from your lift distribution
      real_lift = rlnorm(n(), latent_lift, latent_sd),
      # Compute the power to detect that lift given the sample size you have
      actual_power = map2_dbl(real_lift, n_per_group, one_sided_power),
      # Simulate detecting the lift
      detect = as.logical(rbinom(n(), 1, actual_power)),
      # Did you implement the result or not?
      # If you didn't, this is equivalent to a lift of 1
      # and won't change the product.
      result = if_else(detect, real_lift, 1),
    ) %>% 
    group_by(sim) %>% 
    #finally, take the product, grouping among simulations.
    summarise(lift = prod(result))  %>% 
    ungroup %>% 
    summarise(lift = mean(lift)) %>% 
    pull(lift)
  
})


r <- crossing(
  mde = seq(1.01, 1.1, 0.001),
  n_uuids = 10000000,
  baseline = 0.08,
  latent_lift = 0.01,
  latent_sd = 0.1,
  max_exp_per_year = 12
) %>% 
  mutate(
    n_per_group = samp_size(mde, baseline, n_uuids, max_exp_per_year),
    n_experiments_per_year = floor(n_uuids/(2*n_per_group)),
    lift = simr(n_experiments_per_year, n_per_group, latent_lift, latent_sd)
  ) 


r %>% 
  mutate(n_per_group = log(n_per_group)) %>% 
  select(mde, n_per_group, n_experiments_per_year, lift) %>% 
  pivot_longer(n_per_group:lift) %>% 
  ggplot(aes(mde, value)) + 
  geom_step() + 
  facet_wrap(~name, scales='free', ncol=1)
