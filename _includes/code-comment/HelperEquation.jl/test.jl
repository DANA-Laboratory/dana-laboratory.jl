module test
  #ماژول های مورد استفاده
  using HelperEquation #توابع لازم برای تحلیل مدل
  using IdealGasEos #مدل گاز کامل
  using CpIdeal #ظرفیت حرارتی گاز کامل
  using Calculus #ماژول شامل توابع کلی در جبر
  function testforIdealGasModelWithCp()
    #یک مدل جدید از گاز کامل ساخته شود
    DNIdel=DANAIdealGasEos()
    #دو مقدار از متغییرها مشخص گردد
    DNIdel.P=12.0
    DNIdel.T=120.0
    #پارامتر نوع ماده
    DNIdel.CASNO="95-63-6" #Trimethylbenzene
    #پارامتر نوع تقریب ظرفیت گرمایی
    DNIdel.usePolynomialEstimationOfCp=true
    #از ماژول ظرفیت حرارتی ضرایب برای این گاز استخراج میشود
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = C0Poly("95-63-6")
    #از ماژول مدل گاز کامل معادلات حاکم تعیین می شود
    setEquationFlow(DNIdel)
    #از ماژول تحلیل مدل متغییرها و پارامترهای معلوم جایگذاری میشوند
    a=replace(DNIdel)
    #معادله 1 چاپ شود
    println(a[1]) #=> :(-(*(12.0,v),*(8314.4621,120.0)))
    #کلیه معادلات ساده شود 
    a=map(simplify,a)
    #معادله 1 چاپ شود
    println(a[1]) #=> :(-(*(12.0,v),997735.452))
    #یک دستگاه معادلات جبری تشکیل شود
    vals,vars=analysis(a)
    println(a) #=>
    #:(-(*(12.0,v),997735.452))
    #:(-(Cp,78910.79999999999))
    #:(-(ICpOnTDT,211746.4556136655))
    #:(-(Cv,-(Cp,8314.4621)))
    #:(-(ICpDT,6.785928e6))
    #:(-(u,-(ICpDT,997735.452)))
    #:(-(h,ICpDT))
    #:(-(s,-(ICpOnTDT,20660.662161700304)))
    println(vals) #=>
#12      0       0       0       0       0       0       0       -997735.452
#0       1       0       0       0       0       0       0       -78910.79999999999
#0       0       1       0       0       0       0       0       -211746.4556136655
#0       -1      0       1       0       0       0       0       8314.4621
#0       0       0       0       1       0       0       0       -6785928
#0       0       0       0       -1      1       0       0       997735.452
#0       0       0       0       -1      0       1       0       0
#0       0       -1      0       0       0       0       1       20660.662161700304
    println(vars) #=>
#v
#Cp
#ICpOnTDT
#Cv
#ICpDT
#u
#h
#s
#constant
    #محاسبه ماتریس سطری پلکانی کاهش یافته -گوس جردن
    rVls=rref(vals)
    #ستون آخر ماتریس سطری پلکانی مقدار مجهولات است
    DNIdel.v=-1*last(rVls[1,:])
    DNIdel.Cp=-1*last(rVls[2,:])
    println(rVls) #=>
# 1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  -83144.6
# 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  -78910.8
# 0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  -211746.0
# 0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  -70596.3
# 0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  -6.78593e6
# 0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  -5.78819e6
# 0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  -6.78593e6
# 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  -191086.0
  end
  function testUpdate()
    DNIdel=DANAIdealGasEos()
    DNIdel.P=12.0
    DNIdel.T=120.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = C0Poly("95-63-6")
    setEquationFlow(DNIdel)
    rVls,vars=solve(DNIdel)
    update!(DNIdel,rVls,vars)
    a=replace(DNIdel)
    println(map(eval,a))
  end
  function testIDealGas()
    ######Temprature is undef#######
    DNIdel=DANAIdealGasEos()
    DNIdel.P=2000.0
    DNIdel.v=4000.0
    DNIdel.CASNO="95-63-6"
    DNIdel.usePolynomialEstimationOfCp=true
    DNIdel.C1,DNIdel.C2,DNIdel.C3,DNIdel.C4,DNIdel.C5 = C0Poly("95-63-6")
    setEquationFlow(DNIdel)
    somthingUpdated=true
    fullDetermined=false
    while (somthingUpdated && !fullDetermined)
      rVls,vars=solve(DNIdel)
      somthingUpdated,fullDetermined=update!(DNIdel,rVls,vars)
    end
    dump(DNIdel)
    #a=replace(DNIdel)
    #println(a)
  end
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
end