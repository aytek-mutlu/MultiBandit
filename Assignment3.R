

#Gittins index calculation function
Gittins <- function(alpha,N){
  m = matrix(c(0),N,N)
  for(a in 1:N){
    for(b in 1:N){
      V = matrix(c(0),N,N)
      diff = 1.0
      while(diff > 0.001){
        V_temp  = V
        for(k in a:N){
          for(l in b:N){
            p = (k)/(k+l)
            V[k,l] = max(p + alpha * (p * V_temp[min(k+1,N),l] + (1-p) * V_temp[k,min(l+1,N)]), V_temp[a,b])
          } 
        }
        diff = V[a,b] - V_temp[a,b]
      }
      m[a,b] = V[a,b]*(1-alpha)
    }
  }
return(m)
}

FullInformation <- function(alpha,N){
  p_bandits = runif(2)
  arm = which.max(p_bandits)
  r=0
  for(t in 1:N){
    p_arm = p_bandits[arm]
    x = rbinom(1,1,p_arm)
    r = r + alpha^(t-1) *x
      }
  return(r)
}

#Bayesian with Gittins Index
Bayesian <- function(alpha,N){
  arm_states = matrix(c(0),2,2)
  p_bandits = runif(2)
  arm = which.max(p_bandits)
  r=0
  for(t in 1:N){
    p_arm = (arm_states[arm,2]+1)/(rowSums(arm_states)[arm]+2)
    x = rbinom(1,1,p_arm)
    arm_states[arm,x+1] = arm_states[arm,x+1]+1
    r = r + alpha^(t-1) *x
    arm = ifelse(gittins_table[arm_states[1,2]+1,arm_states[1,1]+1] > gittins_table[arm_states[2,2]+1,arm_states[2,1]+1],1,2)
  }
  return(r)
}

#Thompson Sampling
Thompson <- function(alpha,N){
  arm_states = matrix(c(0),2,2)
  p_bandits = runif(2)
  arm = which.max(p_bandits)
  r=0
  for(t in 1:N){
    p_arm = (arm_states[arm,2]+1)/(rowSums(arm_states)[arm]+2)
    x = rbinom(1,1,p_arm)
    arm_states[arm,x+1] = arm_states[arm,x+1]+1
    r = r + alpha^(t-1) *x
    arm = ifelse(rbeta(n=1,shape1=arm_states[1,2]+1,shape2 =arm_states[1,1]+1)>rbeta(n=1,shape1=arm_states[2,2]+1,shape2 =arm_states[2,1]+1),1,2)
  }
  return(r)
}

#Greedy Q-Learning
Greedy <- function(alpha,N,epsilon,initial_prob){
  
  Q = matrix(c(initial_prob),1,2)
  p_bandits = runif(2)
  r=0
  
  for(t in 1:N){
    
    if(runif(1)<epsilon)
      arm = rbinom(1,1,0.5)+1
    else
      arm = which.max(Q)
    
    x = rbinom(1,1,p_bandits[arm])
    r = r + alpha^(t-1) *x
    
    Q[arm] = (1/t) * alpha^(t-1) * x + (1-(1/t))*Q[arm]
  }
  return(r)
}

#calculate Gittins Index
alpha = 0.95
N = 40
gittins_table = Gittins(alpha,N)

#Simulation
N_sim = 1000
full_info_reward = 0
bayesian_reward = 0
thompson_reward = 0
greedy_reward = 0
eps_greedy_reward = 0
opt_greedy_reward = 0

for(n in 1:N_sim){
  
  fi = FullInformation(alpha,N-1)
  ba = Bayesian(alpha,N-1)
  th = Thompson(alpha,N-1)
  gr = Greedy(alpha,N-1,epsilon=0,initial_prob = 0)
  opt_gr = Greedy(alpha,N-1,epsilon=0,initial_prob = 1)
  eps_gr = Greedy(alpha,N-1,epsilon=0.3,initial_prob = 0)
  

  full_info_reward = full_info_reward + fi
  bayesian_reward = bayesian_reward + ba
  thompson_reward = thompson_reward + th
  greedy_reward = greedy_reward + gr
  opt_greedy_reward = opt_greedy_reward + opt_gr
  eps_greedy_reward = eps_greedy_reward + eps_gr
}

full_info_reward = full_info_reward/N_sim
bayesian_reward = bayesian_reward/N_sim
thompson_reward = thompson_reward/N_sim
greedy_reward = greedy_reward/N_sim
opt_greedy_reward = opt_greedy_reward/N_sim
eps_greedy_reward = eps_greedy_reward/N_sim

print(c("Full info:",full_info_reward))
print(c("Bayesian:",bayesian_reward))
print(c("Thompson Sampling:",thompson_reward))
print(c("Greedy Q-Learning:",greedy_reward))
print(c("Optimistic Greedy Q-Learning:",opt_greedy_reward))
print(c("Epsilon-Greedy Q-Learning:",eps_greedy_reward))

#Find optimal epsilon
epsilons = c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
eps_greedy_rewards = c()

for(eps in epsilons){
  eps_greedy_reward = 0
  for(n in 1:N_sim){
    eps_gr = Greedy(alpha,N-1,epsilon=eps,initial_prob = 0)
    eps_greedy_reward = eps_greedy_reward + eps_gr
  }
  eps_greedy_rewards = c(eps_greedy_rewards,eps_greedy_reward)
}
eps_greedy_rewards = eps_greedy_rewards/N_sim


plot(epsilons,eps_greedy_rewards,type='o',xlab='Epsilons',ylab='Discounted No. of Successes')





