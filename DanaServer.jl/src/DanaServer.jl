module DanaServer
  import DanaTypes.Dana,DanaTypes._Integer,DanaTypes._Switcher,DanaTypes._Real,DanaTypes._Boolean
  using DanaModels
  using JSON
  export createServer

  jsonResponce(data::Dana)=json(data.immute)
  jsonResponce(data::Any)="not defined JsonResponse for DataType: "*string(typeof(data))
  
  function getModuleNames(mod::Module)
    ret=(String=>Any)[]
    mn=names(mod)
    for m in mn
      mm=getfield(mod,m)
      if isa(mm,Module) && !is(mm,mod)
        ret[string(mm)]=getModuleNames(mm)
      elseif (isa(mm,DataType) || isa(mm,Function))
        ret[string(m)]=nothing
      end
    end
    return ret
  end
  function createServer(serverport::Int)
    danaModelsJsonNames=json(getModuleNames(DanaModels))
    c = Base.Condition()
    @async begin  
      try
        server = listen(serverport)
        Base.notify(c)
        println("Listen to ",serverport)
        while true
          sock = accept(server)
          @async while true
            req="DC\r\n"
            try
              req=readline(sock)
            end
            if req=="DC\r\n"
              print("DC request\n")
              break
            elseif req=="DanaModels\r\n"
              print("DanaModels request\n")
              println(sock,danaModelsJsonNames)
            else
              try
                epreq=eval(parse(req))
                if isa(epreq,DataType)
                  # check if requesr is Model or DataType
                  isReqModel::Bool=true
                  if searchindex(req,".EMLtypes.")>0 
                    println("DataType request:",req) 
                    isReqModel=false
                  else
                    println("Model request:",req)
                  end
                  jresp::String=""
                  # inheriented datatype
                  if haskey(Main.NamesOfTypes.inventory,epreq)
                    jresp="{names:"*json(map(string,collect(keys(Main.NamesOfTypes.inventory[epreq]))))*",types:"*json(map(string,collect(values(Main.NamesOfTypes.inventorytyp[epreq]))))*"}"
                    println("responce size=",write(sock,jresp))
                  else # base datatype
                    if isReqModel
                      tmp=epreq()
                      jresp="{names:"*json(map(string,names(tmp)))*"}"
                    else
                      tmp=epreq()
                      jresp=jsonResponce(tmp)
                    end
                    println(write(sock,jresp))
                  end
                else
                  println(sock,"unknown request")
                  print("unknown request: ",req)
                end
              catch 
                println(sock,"error parsing unknown request")
                print("error parsing unknown request: ",req)
              end
            end
          end
        end
      catch
        println("Error")
        Base.notify(c)
      end
    end
    wait(c)
  end
end