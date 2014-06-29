# REF[1] Engineering and Chemical Thermodynamics, 2nd Edition, Milo D. Koretsky
# REF[2] http://en.wikipedia.org/wiki/Departure_function
module PengRobinson
  # Units J,Kmol,Kelvin,pascal
  using DanaTypes
  export DANAPengRobinson,setEquationFlow
  type  DANAPengRobinson <: DanaModel
      DANAPengRobinson()=begin
        new(8314.4621,pi,"",NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
          [
						:(teta=acos(r/q^1.5)),
						:(Z1=-2*sqrt(q)*cos(teta/3)-beta/3),
						:(Z2=-2*sqrt(q)*cos((teta+2*pi)/3)-beta/3),
						:(Z3=-2*sqrt(q)*cos((teta+4*pi)/3)-beta/3),
						:(P=R*T/(v-b)-(d*(1+k*(1-sqrt(T/Tc)))^2)/(v*v+2*b*v-b*b)),
            :(Z=P*v/R/T),
						:(k=0.37464+1.54226*af-0.26992*af^2),
						:(A=(d*(1+k*(1-sqrt(T/Tc)))^2)*P/R^2/T^2),
						:(b=0.0778*R*Tc/Pc),
						:(B=b*P/R/T),
						:(beta=B-1),
						:(gama=A-3*B^2-2*B),
						:(delta=B^3+B^2-A*B),
						:(q=(beta*beta-3*gama)/9),
						:(r=(2*beta^3-9*beta*gama+27*delta)/54),
            :(h_Dep=R*Tc*((T/Tc)*(Z-1)-2.078*(1+k)*sqrt((1+k*(1-sqrt(T/Tc)))^2)*log((Z+2.414*B)/(Z-0.414*B)))), #REF[2]
            :(s_Dep=R*(log(Z-B)-2.078*k*((1+k)/sqrt(T/Tc)-k)*log((Z+2.414*B)/(Z-0.414*B)))), #REF[2]
						:(h_Dep2=((-4*(b^3*R*T*Tc-2*b^2*R*T*Tc*v+d*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))*v^2-b*v*(d*(Tc-2*k*(-1+sqrt(T/Tc))*Tc+k^2*(T+Tc-2*sqrt(T/Tc)*Tc))+R*T*Tc*v)))/(Tc*(b-v)*(b^2-2*b*v-v^2))-(sqrt(2)*d*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(-1+(b+v)/(sqrt(2)*b)))/b+(sqrt(2)*d*(1+k)*(-1+k*(-1+sqrt(T/Tc)))*log(1+(b+v)/(sqrt(2)*b)))/b)/4),
            :(d=0.45724*R^2*Tc^2/Pc),
						:(o=cbrt((r^2-q^3)^0.5+abs(r))),
            :(Z=-sign(r)*(o+q/o)-beta/3)
						#:(PTv=-(((11431*k^2*Tc*v-11431*b*k^2*Tc)*R^2+(-25000*Pc*v^2-50000*b*Pc*v-25000*b^2*Pc)*R)*sqrt(T)+sqrt(Tc)*((-11431*k^2-11431*k)*Tc*v+(11431*b*k^2+11431*b*k)*Tc)*R^2)/((25000*Pc*v^3+25000*b*Pc*v^2-75000*b^2*Pc*v+25000*b^3*Pc)*sqrt(T))),
						#s calculation REF[1] EQ(5.31)
						#:(PTvIv=R*(sqrt((2)*(50000*b*Pc*sqrt(T)-11431*k*R*(-(k*sqrt(T))+sqrt(Tc)+k*sqrt(Tc))*Tc)*atanh((b+v)/(sqrt(2)*b))+25000*b*Pc*sqrt(T)*(4*log(-b+v)-log(-b^2+2*b*v+v^2))))/(50000*b*Pc*sqrt(T))),
						#u calculation REF[1] EQ(5.36)
						#:(TmuPTvmiPIV=-(P*v)+(R*sqrt(T)*(50000*b*Pc*Sqrt(T)+11431*k^2*R*Sqrt(T)*Tc-11431*k*R*Tc^(3/2)-11431*k^2*R*Tc^(3/2))*atanh((b+v)/(sqrt(2)*b)))/(25000*sqrt(2)*b*Pc)+2*R*T*log(-b+v)-(R*T*log(-b^2+2*b*v+v^2))/2),
						#:(PvT=-(((-11431*k^2*Tc*v^3+11431*b*k^2*Tc*v^2+11431*b^2*k^2*Tc*v-11431*b^3*k^2*Tc)*R^2+(12500*Pc*v^4+50000*b*Pc*v^3+100000*b^2*Pc*v^2-62500*b^4*Pc)*R)*T+sqrt(Tc)*((22862*k^2+22862*k)*Tc*v^3+(-22862*b*k^2-22862*b*k)*Tc*v^2+(-22862*b^2*k^2-22862*b^2*k)*Tc*v+(22862*b^3*k^2+22862*b^3*k)*Tc)*R^2*sqrt(T)+((-11431*k^2-22862*k-11431)*Tc^2*v^3+(11431*b*k^2+22862*b*k+11431*b)*Tc^2*v^2+(11431*b^2*k^2+22862*b^2*k+11431*b^2)*Tc^2*v+(-11431*b^3*k^2-22862*b^3*k-11431*b^3)*Tc^2)*R^2)/(12500*Pc*v^6+25000*b*Pc*v^5-62500*b^2*Pc*v^4-50000*b^3*Pc*v^3+137500*b^4*Pc*v^2-75000*b^5*Pc*v+12500*b^6*Pc)),
						#:(PvTIT=(R*(b+v)*(11431*b^2*(1+k)^2*R*Tc^2*T-22862*b*(1+k)^2*R*Tc^2*v*T+11431*(1+k)^2*R*Tc^2*v^2*T-(45724*b^2*k*(1+k)*R*Tc^(3/2)*T^(3/2))/3+(91448*b*k*(1+k)*R*Tc^(3/2)*v*T^(3/2))/3-(45724*k*(1+k)*R*Tc^(3/2)*v^2*T^(3/2))/3+31250*b^3*Pc*T^2+(b^2*(11431*k^2*R*Tc-62500*Pc*v)*T^2)/2+(v^2*(11431*k^2*R*Tc-12500*Pc*v)*T^2)/2-b*v*(11431*k^2*R*Tc+18750*Pc*v)*T^2))/(12500*Pc*(b^3-3*b^2*v+b*v^2+v^3)^2)),
					],Array(Expr,0)
        )
      end
      #paremeters
      R::Float64
      pi::Float64
      CASNO::String
      #variables
      v::Float64
      T::Float64
      Tc::Float64
      P::Float64
      Pc::Float64
			k::Float64
			A::Float64
			b::Float64
			B::Float64
			beta::Float64
			gama::Float64
			delta::Float64
			q::Float64
			r::Float64
			teta::Float64
			af::Float64
			Z1::Float64
			Z2::Float64
			Z3::Float64
			Z::Float64
      o::Float64
      h_Dep::Float64
      s_Dep::Float64
      h_Dep2::Float64
      d::Float64
			#h_Dep_aletr::Float64
      #equations
      equations::Array{Expr,1}
      equationsFlow::Array{Expr,1}
  end
  function setEquationFlow(this::DANAPengRobinson)
    if !isnan(this.q) && !isnan(this.r) && (this.q^3-this.r^2)<0
      this.equationsFlow=this.equations[5:21];
    else
      this.equationsFlow=this.equations[1:19];
    end
    if isnan(this.Z) && !isnan(this.Z1) && !isnan(this.Z2) && !isnan(this.Z3)
      this.Z=max(this.Z1,this.Z2,this.Z3)
    end
  end
end