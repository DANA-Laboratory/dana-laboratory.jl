# REF[1] Chemical Process Design and Integration By Robin Smith
# REF[2] Engineering and Chemical Thermodynamics By Milo D. Koretsky
# REF[3] Perry HandBook
reload ("HelperEquation.jl")
reload ("Tables.jl")
reload ("PengRobinson.jl")
reload ("Calculus.jl")
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testPR()
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testDeparture()
#reload("PengRobinson.jl/test/test_PengRobinson.jl");test_PengRobinson.testVariousKnowns()
module test_PengRobinson
  using HelperEquation
  using PengRobinson
  using Tables
  using Roots
  reload ("IdealGasEos.jl/test/test_IdealGasEos.jl")
  using test_IdealGasEos
	function testVariousKnowns()
		#P & T
		# butane
		PR=DANAPengRobinson()
		PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals","106-97-8")
		PR.P=9.47*1e5
		PR.T=80+273.15
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(PR);
      rVls,vars,nonliFuns,nonliVars=solve(PR)
      somthingUpdated,fullDetermined=update!(PR,rVls,vars)
    end
		if fullDetermined
			println("solved! PR.v=",PR.v)
			v=PR.v
		end
		#v & T
		PR=DANAPengRobinson()
		PR.Tc,PR.Pc,PR.af=getValueForCasNo("Criticals","106-97-8")
		PR.v=v
		PR.T=80+273.15
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(PR);
      rVls,vars,nonliFuns,nonliVars=solve(PR)
      somthingUpdated,fullDetermined=update!(PR,rVls,vars)

			println(rVls)
			println(vars)

		end
		if fullDetermined
			println("solved! PR.P=",PR.P)
		end
	end 
  # Verification REF[2] Example 5.4
  function testDeparture()
    DNpr1=DANAPengRobinson()
		DNpr2=DANAPengRobinson()
		DNpr1.P=9.47*1e5
		DNpr1.T=80+273.15
		DNpr2.P=18.9*1e5
		DNpr2.T=120+273.15
    # butane
		tc,pc,af=getValueForCasNo("Criticals","106-97-8")
		DNpr1.Tc=tc
		DNpr1.Pc=pc
		DNpr1.af=af
		DNpr2.Tc=tc
		DNpr2.Pc=pc
		DNpr2.af=af
		somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(DNpr1);
      rVls,vars,nonliFuns,nonliVars=solve(DNpr1)
      somthingUpdated,fullDetermined=update!(DNpr1,rVls,vars)
    end
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      setEquationFlow(DNpr2);
      rVls,vars,nonliFuns,nonliVars=solve(DNpr2)
      somthingUpdated,fullDetermined=update!(DNpr2,rVls,vars)
    end
    #ideal gas solution 4738.739992314279
    println("result is:",4738.739992314279-DNpr1.h_Dep/1000.0+DNpr2.h_Dep/1000.0,"j/mol REF[2] p285 :3522")
    # REF[2] dh=3.522 (kj/mol)
  end
  function testPR()
    h_Dep::Array{Float64,1}=Array(Float64,0)
    v_Calc::Array{Float64,1}=Array(Float64,0)
		###### verification: check thermodynamic prop of acetone Ref[1]:Table(2-185) for p=0.1Mpa #######
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
      push!(v_Calc,DNpr.v)
      ii=ii+1
    end
    # ideal gas
    #println(h_Dep)
    ideal_h=test_IdealGasEos.forAcetone()
    res1=(ideal_h/1000.0+h_Dep/1000.0)
    println("ideal gas dh (calculated)-> ",round((ideal_h-ideal_h[1])/1000.0)[2:end]," (j/mol)")
    println("real gas dh (claculated)->  ",round(res1-res1[1])[2:end]," (j/mol)")
    println("reference values->          ",round(ref_h-ref_h[1])[2:end]," (j/mol)")
	end
end