mutable struct Strategy
  S
  news::Bool
  sl::Int
  
  function Strategy()
    strat = new()
    
  end
end

function enter_sba(strat::Strategy, atS::Int, x::PBWAlgebraElem, atR::Int)
  strat.news = true
  
  insert!(strat.S, atS, x)
end

