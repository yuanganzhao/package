function [t]=armijo(x,dir,grad,fobj,f0)
dtg=dir'*grad; alpha=0.25; beta = 0.5;
for iter = 0:100,
  t=beta^iter; 
  xn=x+t*dir; 
  fn = feval(fobj,xn);
  if(f0 - fn  >= -alpha * t * dtg ),
  break;
  end
end

function [fobj,grad] = ComputeObj(x,A,b)
fobj = 0.5*x'*A*x+b'*x; grad = A*x+b;
 

