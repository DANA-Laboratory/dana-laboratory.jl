#-------------------------------------------------------------------
#* EMSO Model Library (EML) Copyright (C) 2004 - 2007 ALSOC.
#*
#* This LIBRARY is free software; you can distribute it and/or modify
#* it under the therms of the ALSOC FREE LICENSE as available at
#* http://www.enq.ufrgs.br/alsoc.
#*
#* EMSO Copyright (C) 2004 - 2007 ALSOC, original code
#* from http://www.rps.eng.br Copyright (C) 2002-2004.
#* All rights reserved.
#*
#* EMSO is distributed under the therms of the ALSOC LICENSE as
#* available at http://www.enq.ufrgs.br/alsoc.
#*
#*--------------------------------------------------------------------
#* Models to simulate a power plant.
#*--------------------------------------------------------------------
#* Author: Argimiro R. Secchi
#* $Id: power_plant.mso 195 2007-03-07 20:30:12Z arge $
#*-------------------------------------------------------------------
# Declaracao de tipos
export Dif_Temp
typealias DanaDif_Temp DanaRealParametric
type _Dif_Temp
	function _Dif_Temp(_::Dict{Symbol,Any})
		fields::Dict{Symbol,Any}=(Symbol=>Any)[]
		fields[:Default]=0
		fields[:Lower]=-300
		fields[:Upper]=300
		fields[:Unit]="K"
		drive!(fields,_)
		new(fields)
	end
	value::Dict{Symbol,Any}
end
typealias Dif_Temp DanaRealParametric{_Dif_Temp}