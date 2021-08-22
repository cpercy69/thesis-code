rm(list = ls())
#.libPaths("~/R/library_3.6.2")
#.libPaths()  # confirm set correctly
#devtools::install_github("boettiger-lab/sarsop@0.5.0")
library(sarsop)
library(tidyverse)
library(parallel)
library(gridExtra)
library(tictoc)
library(optparse)

tic()
options(mc.cores=parallel::detectCores())

# DOESN'T LIKE LONG DIRECTORY
log_dir <- "/home/n9427759/thesis-code"

parseArgs <- function(){
  # Parses in command line arguments
  #
  # Args:
  #   NA
  #
  # Returns:
  #   Dataframe with the columns ntrees filled with
  #   value supplied by the command line
  #Load in the arguments from the command line
  option_list = list(
    make_option(c("-g", "--growthrate"), type="double", default=0.75,
                help="intrinsic growth rate"));
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser);
  return(args)
}

args <- parseArgs()

### DEFINITIONS: POPULATION MODEL AND UTILITY FUNCTION ###
r <- args$growthrate
K <- 1

f <- function(x, h){
  s <- pmax(x - h, 0) # compares x-h to 0 and outputs which ever is higher
  s + s * r * (1 - s / K)
}

reward_fn <- function(x, h) pmin(x, h)
discount <- 0.95


### CALCULATING MSY, TAC, AND CE ###

S_star <- optimize(function(x) - f(x,0) + x / discount,
                  c(0, 2*K))$minimum

B_MSY <- f(S_star,0)

MSY <- B_MSY - S_star

F_MSY <- MSY / B_MSY
F_TAC <- 0.8 * F_MSY


MSY_policy <- function(x) F_MSY * x
TAC_policy <- function(x) F_TAC * x

escapement_policy <- function(x) pmax(x - S_star,0)

x0 <- K/6
Tmax = 100
do_det_sim <- function(policy, f, x0, Tmax){
  action <- state <- obs <- as.numeric(rep(NA,Tmax))
  state[1] <- x0
  
  for(t in 1:(Tmax-1)){
    action[t] <- policy(state[t])
    obs[t] <- state[t] - action[t]
    state[t+1] <- f(state[t], action[t])
  }
  data.frame(time = 1:Tmax, state, action, obs)
}

det_sims <- 
  list(MSY = MSY_policy,
       TAC = TAC_policy,
       CE = escapement_policy) %>%
  map_df(do_det_sim,
         f, x0, Tmax,
         .id = "method")

write_csv(det_sims, "/home/n9427759/thesis-code/det_sims.csv")


det_sims %>%
  ggplot(aes(time, state, col=method)) +
  geom_line(lwd=1) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "bottom") +
  ylab("Mean biomass")


### SOLVING THE POMDP SOLUTION ###

## Discretize space
states <- seq(0,2, length=70)
actions <- states
observations <- states

index <- function(x, grid) map_int(x, ~ which.min(abs(.x - grid)))

policies <- data.frame(
  CE = index(escapement_policy(states), actions),
  MSY = index(MSY_policy(states), actions),
  TAC = index(TAC_policy(states), actions))


### DEFINE THE POMDP MATRICES ###

meta <- expand.grid(sigma_g = c(0.02, 0.1, 0.15),
                    sigma_m = c(0, 0.1, 0.15),
                    stringsAsFactors = FALSE) %>%
        mutate(scenario = as.character(1:length(sigma_m)))

meta <- meta[9,]
write_csv(meta, file.path(log_dir, "meta.csv"))

models <- parallel::mclapply(1:dim(meta)[1],
                    function(i){
          fisheries_matrices(
          states = states,
          actions = actions,
          observed_states = observations,
          reward_fn = reward_fn,
          f = f,
          sigma_g = meta[i,"sigma_g"][[1]],
          sigma_m = meta[i,"sigma_m"][[1]],
          noise = "normal")
          })

### COMPUTE THE POMDP SOLUTION ###

dir.create(log_dir, FALSE)
system.time(
  alphas <- 
    parallel::mclapply(1:length(models),
    function(i){
      log_data <- data.frame(model = "gs",
                             r = r,
                             K = K,
                             sigma_g = meta[i,"sigma_g"][[1]],
                             sigma_m = meta[i,"sigma_m"][[1]],
                             noise = "normal",
                             scenario = meta[i, "scenario"][[1]])
      
      sarsop(models[[i]]$transition,
             models[[i]]$observation,
             models[[i]]$reward,
             discount = discount,
             precision = 1e-3,
             timeout = 150,
             log_dir = log_dir,
             log_data = log_data)
    })
)


### SIMULATIONS ###

## Simulating the static policies under uncertainty
set.seed(12345)

Tmax <- 100
x0 <- which.min(abs(K/6 - states))
reps <- 100
static_sims <-
  map_dfr(models, function(m){
    do_sim <- function(policy) sim_pomdp(
              m$transition, m$observation, m$reward, discount,
              x0 = x0, Tmax = Tmax, policy = policy, reps = reps)$df
    map_dfr(policies, do_sim, .id = "method")
  }, .id = "scenario")

## Simulating the POMDP policies under uncertainty
set.seed(12345)

unif_prior <- rep(1, length(states)) / length(states)
pomdp_sims <-
  map2_dfr(models, alphas, function(.x, .y){
            sim_pomdp(.x$transition, .x$observation, .x$reward, discount,
                      unif_prior, x0 = x0, Tmax = Tmax, alpha = .y,
                      reps = reps)$df %>%
            mutate(method = "POMDP")
          },
          .id = "scenario")

write_csv(pomdp_sims, file.path(log_dir, "pomdp_sims.csv"))

pomdp_sims <- read_csv(file.path(log_dir,"pomdp_sims.csv")) %>%
  mutate(scenario = as.character(scenario))

sims <- bind_rows(static_sims, pomdp_sims) %>%
  left_join(meta) %>%
  mutate(state = states[state], action = actions[action]) %>%
  select(time, state, rep, method, sigma_m, sigma_g, value)

write_csv(sims, file.path(log_dir, "sims.csv"))

## Figure S2: Full Simulation Plots
sims <- read.csv(file.path(log_dir,"sims.csv"), col_types = "inicnnn")

sims %>%
  select(time, state, rep, method, sigma_m, sigma_g) %>%
  group_by(time, method, sigma_m, sigma_g) %>%
  summarise(mean = mean(state), sd = sd(state)) %>%
  ggplot(aes(time, mean, col=method, fill=method)) +
  geom_line() +
  geom_ribbon(aes(ymax = mean + sd, ymin = mean-sd), col = NA, alpha = 0.1) +
  facet_grid(sigma_m ~ sigma_g, labeller = label_bquote(sigma[m] == 0.25,
                                                        sigma[g] == 0.02)) + 
  coord_cartesian(ylim = c(0, 1))




















