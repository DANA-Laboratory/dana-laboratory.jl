reload ("HelperEquation.jl")
reload ("Tables.jl")
reload ("PengRobinson.jl")
reload ("Calculus.jl")
module test
  using HelperEquation
  using PengRobinson
  using Tables
  using Roots
  function testPR()
		###### verification: check thermodynamic prop of acetone Ref[1]:Table(2-185) #######
    ref_valuse=[25.930,28.002,32.563,36.923,41.200,45.437]
    ii=1
		for T in [328.84,350,400,450,500,550]
			DNpr=DANAPengRobinson();
			DNpr.T=T;
			DNpr.P=0.1e6;
			# acetone
			tc,pc,af=getValueForCasNo("Criticals","67-64-1");
			DNpr.Tc=tc;
			DNpr.Pc=pc;
			DNpr.af=af;
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
      println("for temperature=",DNpr.T," volume/kmol=",DNpr.v," reference_value=",ref_valuse[ii]) #Table[2.185]
      ii=ii+1
    end
	end
end