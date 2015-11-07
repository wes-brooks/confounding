tol <- 0.1
diff <- 1
f.old <- obj$fn(obj$par)
f.new <- Inf
par <- obj$par
t <- 1e-2

while(diff>tol) {
  while(f.new > f.old) {
    gr.new <- obj$gr(par)
    
    if (obj$fn(par + t*gr.new))
  }
}