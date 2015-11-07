conjugate.gradient <- function(objective, gradient, y, eta, S, wt, ltau, M, V, verbose=TRUE, tol=1e-4) {
  #Initial parameters:
  n = nrow(S)
  r <- ncol(S)
  
  finished <- FALSE
  f.new <- objective(y, eta, S, wt, ltau, M, V)
  f.old <- Inf
  check <- Inf
  
  while(!finished) {
    
    #These iterations restart conjugacy:
    converged = FALSE
    iter <- 0
    while (!converged && iter<100) {
      iter <- iter+1
      
      #Prepare to iterate conjugate gradient descent:
      f.outer <- f.old
      f.old <- Inf
      t <- 1
      conv.inner <- FALSE
      i <- 0
      
      if (verbose) cat("Iterating to estimate V")
      while(f.new < f.old && !converged && !conv.inner && i<(r*(r+1)/2 - 1)) {
        i <- i+1
        
        dir.new <- gradient(y, eta, S, wt, ltau, V)
        dir.new <- dir.new / sqrt(sum(dir.new^2))
        
        #First iteration, ignore conjugacy - thereafter, use it.
        #s.new is the vector of the new step (in parameter space)
        if (i==1) {  
          s.new <- dir.new
        } else {
          conj <- (sum(dir.new^2) + sum(dir.new * dir.old)) / sum(dir.old^2)
          s.new <- dir.new + conj * s.old
        }
        
        V.step <- matrix(s.new, r, r)
        
        #Find the optimal step size
        #Backtracking: stop when the loss function is majorized
        f.proposed <- objective(y, eta, S, wt, ltau, M, V + t*V.step)
        condition <- (f.proposed > f.new - sum((t*s.new)*dir.new) - 1/(2*t)*sum((t*s.new)^2))[1]
        while(condition && t > .Machine$double.eps) {
          t = 0.5*t
          f.proposed <- objective(y, eta, S, wt, ltau, M, V + t*V.step)
          condition = (f.proposed > f.new - sum((t*s.new)*dir.new) - 1/(2*t)*sum((t*s.new)^2))[1]
          
          #This is the final stopping rule: t gets so small that 1/(2*t) is Inf
          if (is.na(condition)) {
            converged = TRUE
            condition = FALSE
          }
        }
        
        #Find the optimal step
        V.proposed <- V + t*V.step
        
        #Make t a little bigger so next iteration has option to make larger step:
        t = t / 0.5 / 0.5
        
        #save for next iteration:
        dir.old <- dir.new
        s.old <- s.new
        
        #Only save the new parameters if they've decreased the loss function
        if (f.proposed < f.old) {
          V <- V.proposed
        }
        f.old <- f.new
        f.new <- f.proposed
        
        if (verbose) cat(".")
        if ((f.old - f.new) < tol * (tol+abs(f.old))) conv.inner = TRUE
      }
      
      if (verbose) cat("done!\n")
      if (verbose) cat(paste("Iteration: ", iter, "; Objective: ", f.new, "; Inner iterations: ", i, "\n", sep=""))
      if ((f.outer - f.new) < tol * (tol+abs(f.outer))) converged = TRUE
    }

    finished <- TRUE
  }
  
  list(V=V, likelihood=f.new)
}