module IdealGasEos
  #Units J,Kmol,Kelvin
  #انواع داده های داخلی دانا را استفاده میکنیم
  using DanaTypes
  #مدل صادره شامل یک ساختمان داده و رویه های مربوطه میباشد
  export DANAIdealGasEos,setEquationFlow
  #ساختمان داده مدل گاز ایدآل
  type  DANAIdealGasEos <: DanaModel
  #مولد ساختمان
      DANAIdealGasEos()=begin
        #مقادیر عناصر مدل به ترتیب تعریف شدن
        new(8314.4621,"",true,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
          #معادلات مدل
          #معادلات مدل جمعا 11 معادله میباشد که 8 معادله مستقلند
          [
            :(P*v=R*T),#معادله حالت
            :(Cp=C1+C2*T+C3*T^2+C4*T^3+C5*T^4),#Poly Cp
            :(Cp=C1+C2*((C3/T)/sinh(C3/T))^2+C4*((C5/T)/cosh(C5/T))^2),#Hyper Cp
            :(ICpOnTDT=C2*T+(C3*T^2)/2+(C4*T^3)/3+(C5*T^4)/4+C1*log(T)),#Integral of Cp/T Poly
            :(ICpOnTDT=(C2*C3*coth(C3/T)+C1*T*log(T)+C4*T*log(cosh(C5/T))-C2*T*log(Sinh(C3/T))-C4*C5*tanh(C5/T))/T),#Integral of Cp/T Hyper
            :(Cv=Cp-R),#Cv Def
            :(ICpDT=C1*T+1/60*T^2*(30*C2+T*(20*C3+3*T*(5*C4+4*C5*T)))),#Integ of Cp Poly
            :(ICpDT=C1*T+C2*C3*coth(C3/T)-C4*C5*tanh(C5/T)),#Integ of Cp Hyper
            :(u=ICpDT-R*T), #Internal energy dep
            :(h=ICpDT), #Enthalpy def
            :(s=ICpOnTDT-R*log(P)) #Entropy
          ],Array(Expr,0)
        )
      end
      ##################### پارامتر های مدل قبل از حل مشخص میباشند ##################################
      R::Float64 #پارامتر ثابت R
      CASNO::String #کد اختصاری ماده
      usePolynomialEstimationOfCp::Bool #کدام تقریب برای محاسبه ظرفیت گرمایی استفاده شود
      C1::Float64 #ضریب
      C2::Float64 #ضریب
      C3::Float64 #ضریب
      C4::Float64 #ضریب
      C5::Float64 #ضریب
      ##################### متغییر های مدل ##################################
      #تعداد 8
      v::Float64 #حجم مخصوص
      T::Float64 #دما
      P::Float64 #فشار
      Cp::Float64 #ظرفیت گرمایی
      Cv::Float64 #ظرفیت گرمایی حجم ثابت
      u::Float64 #انرژی داخلی
      h::Float64 #آنتالپی
      s::Float64 #آنتروپی
      ##################### متغییرهای وابسته ##################################
      #تعداد 2
      ICpOnTDT::Float64 #Integral of Cp/T Poly
      ICpDT::Float64 #Integral of Cp/T Poly
      #تعداد کل متغییرها 10 میباشد
      equations::Array{Expr,1} #آرایه ای برای نگه داری کلیه معادلات لازم
      equationsFlow::Array{Expr,1} #آرایه ای از معادلات حاکم
  end
  #با توجه به شرایط مدل 8 معادله مستقل از 11 ملادله برگزیده میشود
  function setEquationFlow(this::DANAIdealGasEos)
    if this.usePolynomialEstimationOfCp
      this.equationsFlow=[this.equations[1],this.equations[2],this.equations[4],this.equations[6],this.equations[7],this.equations[9],this.equations[10],this.equations[11]]
    else
      this.equationsFlow=[this.equations[1],this.equations[3],this.equations[5],this.equations[6],this.equations[8],this.equations[9],this.equations[10],this.equations[11]]
    end
  end
end
