#استخراج ضرایب مناسب برای دو روش محاسبه ظرفیت گرمایی با استفاده از تقریب چند جمله ای یا هایپربولیک
module CpIdeal
#توابع صادره از ماژول 
  export C0Poly,C0Hyper
#فایل های جداول باز شوند
  cp_polynomial=open("./share/julia/site/v0.2/CpIdeal.jl/src/Tables/perryHeatCapIdealGas(Table2-155).table");
  cp_hyperbolic=open("./share/julia/site/v0.2/CpIdeal.jl/src/Tables/perryHeatCapIdealGas(Table2-156).table");
#جدول 2-155 هندبوک پری در ماتریس بارگذاری شود
  data_poly,header_poly=readdlm(cp_polynomial,';',has_header=true);
#جدول 2-156 هندبوک پری در ماتریس بارگذاری شود
  data_hyper,header_hyper=readdlm(cp_hyperbolic,';',has_header=true);
#نیاز به فایل نیست
  close(cp_polynomial)
  close(cp_hyperbolic)
#ثابت های تقریف چند جمله ای را برای یک ماده با کد مشخص برمیگرداند
  function C0Poly(CasNo::String)
    length=size(data_poly)[1]
    i=1
    while (i<=length && data_poly[i,4]!=CasNo) 
      i+=1;
    end
    if (i<=length)
      return (data_poly[i,6],data_poly[i,7],data_poly[i,8],data_poly[i,13],data_poly[i,14])
    end
  end
#ثابتهای تقریب هایپربولیک را برای یک ماده با کد مشخص برمیگرداند
  function C0Hyper(CasNo::String)
    length=size(data_hyper)[1]
    i=1
    while (i<=length && data_hyper[i,4]!=CasNo)
      i+=1;
    end
    if (i<=length)
      return (data_hyper[i,6],data_hyper[i,7],data_hyper[i,8],data_hyper[i,9],data_hyper[i,10])
    end
  end
end