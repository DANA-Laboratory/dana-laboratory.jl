reload ("HelperEquation.jl")
reload ("Tables.jl")
reload ("PengRobinson.jl")
reload ("Calculus.jl")
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testPR()
module test_PengRobinson
  using HelperEquation
  using PengRobinson
  using Tables
  using Roots
  reload ("IdealGasEos.jl/test/test_IdealGasEos.jl")
  using test_IdealGasEos
  function testPR()
    h_Dep::Array{Float64,1}=Array(Float64,0)
		###### verification: check thermodynamic prop of acetone Ref[1]:Table(2-185) #######
    ref_h=[47.730,49.643,54.255,59.166,64.436,70.066]*1000
    ref_s=[0.16988,0.17552,0.18783,0.19939,0.21049,0.22122]*1000
    ref_valuse=[25.930,28.002,32.563,36.923,41.200,45.437]
    ii=1
		for T in [328.84,350,400,450,500,550]
			DNpr=DANAPengRobinson()
			DNpr.T=T
			DNpr.P=0.1e6
			# acetone
			tc,pc,af=getValueForCasNo("Criticals","67-64-1");
			DNpr.Tc=tc
			DNpr.Pc=pc
			DNpr.af=af
			somthingUpdated=true
			fullDetermined=false
			nonliFuns::Array{Function,1}=Array(Function,0)
			nonliVars::Array{Array{String,1},1}=Array(Array{String,1},0)
			while (somthingUpdated && !fullDetermined)
				while (somthingUpdated && !fullDetermined)
          setEquationFlow(DNpr);
					rVls,vars,nonliFuns,nonliVars=solve(DNpr)
					somthingUpdated,fullDetermined=update!(DNpr,rVls,vars)
				end
				if !fullDetermined
					i=1
					fullDetermined=true
					while (i<=length(nonliFuns))
						if length(nonliVars[i])==1
							result=Roots.fzero(nonliFuns[i],[0,typemax(Int64)])
							HelperEquation.setfield(DNpr,nonliVars[i][1],result)
							somthingUpdated=true
						else
							fullDetermined=false
						end 
						i=i+1
					end
				end
			end
      #println("T=",DNpr.T," v=",DNpr.v," ref_val=",ref_valuse[ii]) #Table[2.185]
      #println(" Dh=",DNpr.h_Dep," Ds=",DNpr.s_Dep) #Table[2.185]
      push!(h_Dep,DNpr.h_Dep)
      ii=ii+1
    end
    # ideal gas
    println(h_Dep)
    println(test_IdealGasEos.forAcetone())
    println(ref_h-ref_h[1])
	end
end