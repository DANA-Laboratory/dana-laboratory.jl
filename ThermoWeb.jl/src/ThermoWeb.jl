# REF[1] Engineering and Chemical Thermodynamics, 2nd Edition, Milo D. Koretsky
# includes 8 equations and 12 variables
module ThermoWeb
  using DanaTypes
  export DANAThermoWeb,setEquationFlow
  type  DANAThermoWeb <: DanaModel
      DANAThermoWeb()=begin
        new(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
          [ # REF[1] p274
            :(nPvs=Psv*svP),#Cyclic_1 
            :(svP=PTs),#Maxwell_1
            :(PTs=nPsT*sTP),#Cyclic_2
						:(nPsT=TvP),#Maxwell_2
            :(nPvT=TvP*PTv),#Cyclic_3
						:(PTv=svT),#Maxwell_3
						:(nTvs=Tsv*svT),#Cyclic_4
            :(nTvs=Psv) #Maxwell_4
          ],Array(Expr,0)
        )
      end
      #paremeters
      #variables
      nPvs::Float64
      Psv::Float64
      svP::Float64
      PTs::Float64
      nPsT::Float64
      sTP::Float64
      TvP::Float64
			nPvT::Float64
      PTv::Float64
      svT::Float64
      Tsv::Float64
			nTvs::Float64
      #equations
      equations::Array{Expr,1}
      equationsFlow::Array{Expr,1}
  end
  function setEquationFlow(this::DANAIdealGasEos)
		this.equationsFlow=this.equations
  end
end
