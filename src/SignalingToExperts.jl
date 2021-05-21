module SignalingToExperts

export errmktclear, equil, f, indif, mktclear2

#This file create the numerical simulations used to generate figure 6

"""
    We will use the following packages for our replication study


"""
using Roots, Plots, LaTeXStrings, Statistics, ColorSchemes, Calculus, QuadGK

"""
    We will first start with the functions needed to plot the graphs for figure 6. 
        These Functions are: errmktclear, equil, f, indif, mktclear2. 


"""

function errmktclear(pi::Any,theta::Any)


    g(t)= Delta./((lambda-t) * lambdahat/lambda + pi *(1-lambdahat)) *f(t)
   
    y= quadgk(g,theta,lambda) 
      return y[1]-1 
end


function equil(pi::Any)

    theta = indif(pi)

    y = errmktclear(pi,theta)

    return y
    end

function f(theta::Any)

   
    y = A .*(theta>=0) .*(theta<=lambda) .* ((sin((B .*(lambda .- theta) .^C)+E .*pi)+1).^D + F) 
return y
end

function indif(pi::Any) #only one working for now
        
    y = lambda - pi .*lambdahat/lambda .*(1-lambdahat) .*1/(cL/cH-1)   
        
    return y
end

function mktclear2(pi::Any,guesstheta::Any)

    if errmktclear(pi,0) < 0
        y = 0
    else() 
        #guesstheta
        y=find_zero( theta -> errmktclear(pi,theta),(0,guesstheta))
    end 
    
    return y
end 

## Testing these functions gives us the same values as found in the matlab version.
"""
   Now that the functions are generated , we will define the parameters of interest. 

"""
global qH, qL, cH, cL, A, B, C, D, E, F, lambda, lambdahat, Delta  

# Parameters
qH = 1
qL = 0.4

cH = .09
cL = .15

lambda = 0.6
lambdahat = lambda
Delta = 1

##Parameters of F
A = .5
B = 13.3
C = .4
D = 2.8
E = 2.2
F = 0


estar = (qH - qL)/cL
wstar = qH - cH*estar


N = 100
pis = LinRange(0.00001,.99999,N)

N2 = 1000

# Plot density
thetas= zeros(N2)
thetas = LinRange(0,lambda,N2)
plot(thetas,f(thetas))
xlabel("\theta")
ylabel("f[\theta]")

## Plot equilibrium conditions

guesstheta = lambda-0.03

for j = 1:N
    inds[j] = indif(pis[j])
    mcs[j] = mktclear2(pis[j],guesstheta)
    if j>1
        guesstheta = mcs[j-1]
    end
end

plot(pis,inds,"b',pis,mcs,'r','LineWidth",1.2) ## this is all I want to plot

leg = legend("Indifference','Market Clearing")
set(leg,"Interpreter','Latex','Fontsize',15,'Location','NorthEast")
xlabel("$\pi$','Interpreter', 'Latex','FontSize", 20) ## not important
ylabel("$\theta$','Interpreter', 'Latex','FontSize", 20) ## not important


""" Find equilibria near hand-picked points"""
pistar[1] = fzero(equil(pi),0.1) 
pistar[2] = fzero(equil(pi),0.9)

for j = 1:2
    thetastar[j] = mktclear2(pistar[j],guesstheta)
    guesstheta = thetastar[j]
end

## Compute profits in equilibrium & deviation

thetasdev[1,:]=range(thetastar[1],lambda,length=N)
thetasdev[2,:]=range(thetastar[2],lambda,length=N)

profits[1,:] = ((lambda - thetasdev[1,:])*qL + pistar[1]*(1-lambda)*qH)./ (lambda - thetasdev[1,:] + pistar[1]*(1-lambda)) .-wstar

profits[2,:] = ((lambda - thetasdev[2,:])*qL + pistar[2]*(1-lambda)*qH)./ (lambda - thetasdev[2,:] + pistar[2]*(1-lambda)) .-wstar

for j = 1:N

    a= quadgk[ 1 ./ (pistar[1]*(1-lambda)+lambda-t) .*f(t),thetastar[1],thetasdev[1,j]]
    b= quadgk[ 1 ./(pistar[2]*(1-lambda)+lambda-t) .*f(t),thetastar[2],thetasdev[2,j]]

    wD[1,j] = qL + (qH - qL) * (1-cH/cL * a[1]) 

    devprofits[1,j] = qH -wD[1,j]

    wD[2,j] = qL + (qH - qL) * (1-cH/cL * b[1]) 
    
    devprofits[2,j] = qH -wD[2,j]
end

plot(thetasdev[1,:],profits[1,:],thetasdev[1,:],devprofits[1,:])
legend("eq','dev")

plot(thetasdev[2,:],profits[2,:],thetasdev[2,:],devprofits[2,:])
legend("eq','dev")


end
