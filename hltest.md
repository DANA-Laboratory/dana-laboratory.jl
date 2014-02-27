---
layout: default
title: Something with codes
filename: hltest.md
---

```julia
  function testIDealGasForNonlinearSolver()
    ######Temprature is undef#######
    DNIdel=DANAIdealGasEos()
    DNIdel.P=2000.0
    DNIdel.Cp=629657.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = C0Poly("95-63-6")
    setEquationFlow(DNIdel)
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      rVls,vars=solve(DNIdel)
      println("************one solution done************")
      somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
    end
    dump(DNIdel)
    #a=replace(DNIdel)
    #println(a)
  end
```