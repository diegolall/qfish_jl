using LinearAlgebra

function qfish(rho,J,Threshold=1e-14)
    if(minimum(size(rho))==1 || ndims(rho)<2)
        rho = kron(rho,rho')
    end #if 

    rho=rho/tr(rho)

    #make it Hermitian to avoid diag problems
    rho=Hermitian(rho)

    sx,sy=size(rho)
    evals, evecs = LinearAlgebra.eigen(rho)

    if sum(evals .> Threshold) == 1
        #pure state
        #println("Pure state")
        f=tr(J*J*rho)-tr(J*rho)^2
        if abs(imag(f)) < Threshold
            return 4*real(f)
        else
            error("Imaginary part of f is greater than the threshold.")
        end #if
    elseif minimum(evals) > Threshold
        #full rank
        #println("Full rank")
        f=0
        for n in 1:sx
            lambdan=evals[n]
            v=evecs[:,n]'*J
            for m in 1:n-1
                lambdam=evals[m]
                f+=(lambdan-lambdam)^2/(lambdan+lambdam)*abs(v*evecs[:,m])^2
            end #for
        end #for
        if abs(imag(f)) < Threshold
            return 4*real(f)
        else
            error("Imaginary part of f is greater than the threshold.")
        end #if
    else
        #not full rank
        #println("Not full rank")
        f=0
        for n in 1:sx
            lambdan=evals[n]
            v=evecs[:,n]'*J
            for m in 1:n-1
                lambdam=evals[m]
                if abs(lambdan+lambdam) >  2*Threshold
                    f+=(lambdan-lambdam)^2/(lambdan+lambdam)*abs(v*evecs[:,m])^2
                end #if
            end #for
        end #for
        if abs(imag(f)) < Threshold
            return 4*real(f)
        else
            error("Imaginary part of f is greater than the threshold.")
        end #if
    end #if
    error("Unexpected condition: the code should not reach this point in qfish.")
end #function