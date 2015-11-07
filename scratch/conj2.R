#Initial parameters:
n = nrow(data)
p = ncol(data)

finished = FALSE

par <- obj$par
log.lik <- obj$fn
score <- obj$gr

finished <- FALSE
f.new = log.lik(par)
f.old = Inf
tol = sqrt(.Machine$double.eps)
tol = 1e-3
check=Inf

while(!finished) {

  #These iterations restart conjugacy:
  converged = FALSE
  iter <- 0
  while (!converged && iter<100) {
    iter <- iter+1
    
    #Prepare to iterate conjugate gradient descent:
    i <- 0
    f.outer <- f.old
    f.old <- Inf
    t <- 1
    conv.inner <- FALSE
    while(f.new < f.old && !converged && !conv.inner && i<(length(par)-1)) {
      i <- i+1
      
      dir.new <- -score(par)
      dir.new <- dir.new / sqrt(sum(dir.new^2))
      
      #First iteration, ignore conjugacy - thereafter, use it.
      #s.new is the vector of the new step (in parameter space)
      if (i==1) {  
        s.new <- dir.new
      } else {
        conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
        s.new <- dir.new + conj * s.old
      }
      
      #Find the optimal step size
      #Backtracking: stop when the loss function is majorized
      condition = (log.lik(par + t*s.new) > f.new - sum((t*s.new)*dir.new) - 1/(2*t)*sum((t*s.new)^2))[1]
      while(condition && t > .Machine$double.eps) {
        condition = (log.lik(par + t*s.new) > f.new - sum((t*s.new)*dir.new) - 1/(2*t)*sum((t*s.new)^2))[1]
        
        #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
        if (is.na(condition)) {
          converged = TRUE
          condition = FALSE
        }
        
        t = 0.8*t
      }
      
      #Find the optimal step
      step = s.new * t
      
      #Make t a little bigger so next iteration has option to make larger step:
      t = t / 0.8 / 0.8
      p = par + step
      
      #save for next iteration:
      dir.old <- dir.new
      s.old <- s.new
      
      #Only save the new parameters if they've decreased the loss function
      f.proposed <- log.lik(p)
      if (f.proposed < f.old)
        par <- p
      f.old <- f.new
      f.new <- f.proposed
      
      if ((f.old - f.new) < tol * (tol+abs(f.old))) conv.inner = TRUE
    }
    
    if (verbose) cat(paste("Iteration: ", iter, "; Likelihood objective: ", f.new, "; Inner iterations: ", i, "\n", sep=""))
    
    
    if ((f.outer - f.new) < tol * (tol+abs(f.outer))) converged = TRUE
  }
  
  if (tol > .Machine$double.eps) tol = tol/10
  else finished <- TRUE
}