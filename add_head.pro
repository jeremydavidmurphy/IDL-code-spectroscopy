; This routine is used to add a single line to a fits header. As this
; will often be done to multiple fits files (i.e. the same addition to
; many files), it accepts a list of fits files. The string that gets
; added is internal to the code.

; WARNING! This code will overwrite your file with the same fits file
; but new header. So, make sure you're ready before you run this...

;example = 'RA      =  10:59:44.60         / Right ascension'
;example = 'DEC     =  -18:18:33.3         / Declination'
;example = 'EQUINOX =              2000.00 / Equinox of coordinate system'
;example = 'AIRMASS =                 1.55 / Airmass'

PRO add_head, datalist

readcol,datalist,f='a,a',list,galaxy
n0 = n_elements(list)

;*********************************************************************
;addH1 = 'RA      =         04:31:39.80  / Right ascension' ;N1600
;addH2 = 'DEC     =         -05:05:10.0  / Declination' ;N1600
;addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
;addH4 = 'AIRMASS =                1.35  / Airmass'
;addH1 = 'RA      =         09:19:46.80  / Right ascension' ;N2832
;addH2 = 'DEC     =         +33:44:59.0  / Declination' ;N2832
;addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
;addH4 = 'AIRMASS =                1.10  / Airmass'
;addH1 = 'RA      =         12:30:49.40  / Right ascension' ;N4486
;addH2 = 'DEC     =         +12:23:28.0  / Declination' ;N4486
;addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
;addH4 = 'AIRMASS =                1.15  / Airmass'
;addH1 = 'RA      =         05:05:30.60  / Right ascension' ;G191B2B
;addH2 = 'DEC     =         +27:20:17.8  / Declination'  ;G191B2B
;addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
;addH4 = 'AIRMASS =                1.05  / Airmass'
;addH1 = 'RA      =         10:39:36.70  / Right ascension' ;Feige 34
;addH2 = 'DEC     =         +42:06:09.5  / Declination'  ;Feige 34
;addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
;addH4 = 'AIRMASS =                1.05  / Airmass'
addH1 = 'RA      =         07:00:00.00  / Right ascension' ;"standard-feb08"
addH2 = 'DEC     =         +42:00:00.0  / Declination'
addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
addH4 = 'AIRMASS =                1.10  / Airmass'
nstring = 4 ;set this number to the number of lines you're adding (limit 5)
;*********************************************************************

;readcol,datalist,f='a',list

for j=0,n0-1 do begin
    file = list[j]+'.fits'
    data = readfits(file,h)
    nh = n_elements(h)
    if galaxy[j] eq 'n1600' then begin
       addH1 = 'RA      =         04:31:39.80  / Right ascension' ;N1600
       addH2 = 'DEC     =         -05:05:10.0  / Declination'     ;N1600
       addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
       addH4 = 'AIRMASS =                1.35  / Airmass'
    endif
    if galaxy[j] eq 'n2832' then begin
       addH1 = 'RA      =         09:19:46.80  / Right ascension' ;N2832
       addH2 = 'DEC     =         +33:44:59.0  / Declination'     ;N2832
       addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
       addH4 = 'AIRMASS =                1.10  / Airmass'
    endif
    if galaxy[j] eq 'n4486' then begin
       addH1 = 'RA      =         12:30:49.40  / Right ascension' ;N4486
       addH2 = 'DEC     =         +12:23:28.0  / Declination'     ;N4486
       addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
       addH4 = 'AIRMASS =                1.15  / Airmass'
    endif
    if galaxy[j] eq 'g191b2b' then begin
       addH1 = 'RA      =         05:05:30.60  / Right ascension' ;G191B2B
       addH2 = 'DEC     =         +27:20:17.8  / Declination'     ;G191B2B
       addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
       addH4 = 'AIRMASS =                1.05  / Airmass'
    endif
    if galaxy[j] eq 'feige34' then begin
       addH1 = 'RA      =         10:39:36.70  / Right ascension' ;Feige 34
       addH2 = 'DEC     =         +42:06:09.5  / Declination'     ;Feige 34
       addH3 = 'EQUINOX =             2000.00  / Equinox of coordinate system'
       addH4 = 'AIRMASS =                1.05  / Airmass'
    endif
    
    if (nstring eq 1) then hout = [h[0:nh-2],addH1,h[nh-1]]
    if (nstring eq 2) then hout = [h[0:nh-2],addH1,addH2,h[nh-1]]
    if (nstring eq 3) then hout = [h[0:nh-2],addH1,addH2,addH3,h[nh-1]]
    if (nstring eq 4) then hout = [h[0:nh-2],addH1,addH2,addH3,addH4,h[nh-1]]
    if (nstring eq 5) then hout = [h[0:nh-2],addH1,addH2,addH3,addH4,addH5,h[nh-1]]
    if (nstring gt 5) then print,'Yo! Change your code! Too many headers!'
    writefits,file,data,hout
endfor


stop
END

