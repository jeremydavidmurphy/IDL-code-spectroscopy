;----------------------------------------------------------------------
function sign, a
; Implements sign(x)=x/abs(x)
; Michele Cappellari, Leiden, 29 May 2003
return, (a gt 0)*2 - 1
end
;----------------------------------------------------------------------
