reload ("HelperEquation.jl")
reload ("Tables.jl")
module test
  using HelperEquation
  using IdealGasEos
  using Tables
  using Roots
  function testIDealGas()
    ###### verification: check monoxide Enthalpies with Ref[1]:Table(2-180) #######
    DNIdel=DANAIdealGasEos()
    DNIdel.T=298.15
		# monoxide
    DNIdel.CASNO="630-08-0"
    DNIdel.usePolynomialEstimationOfCp=false
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Hyper",DNIdel.CASNO)
    setEquationFlow(DNIdel)
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
			rVls,vars=solve(DNIdel)
      somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
    end
		hst=DNIdel.h;
		ust=DNIdel.u;
		refDelta_h=[-2858,-1692,-1110,-529,54,638,1221,1805,2389,2975,3563,4153,4643,5335,5931,7428,8942,10477,12023,13592,15177,16781,18401,20031,21690,25035,28430,31868,35343,38850,42385,45945,49526,53126,56744,60376,64021,67683,71324,74985,78673,82369,86074,89786,93504,112185,130989,149895,168890];
		i=1;
		for T in [200,240,260,280,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3500,4000,4500,5000]	
			DNIdel=DANAIdealGasEos()
			DNIdel.T=T
			# monoxide
			DNIdel.CASNO="630-08-0"
			DNIdel.usePolynomialEstimationOfCp=false
			DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = getValueForCasNo("C0Hyper",DNIdel.CASNO)
			setEquationFlow(DNIdel)
			somthingUpdated=true
			fullDetermined=false
			while (somthingUpdated && !fullDetermined)
				rVls,vars=solve(DNIdel)
				somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
			end
			println("T=",T," Dh=",(DNIdel.h-hst)/1000," ref value=",refDelta_h[i]," diff=",(DNIdel.h-hst)/1000-refDelta_h[i]);
			i+=1;
		end
	end
end