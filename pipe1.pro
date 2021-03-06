;This is the first piece of the pipeline, which ties together
;radius.pro and radplot.pro, and works from lists. Each pointing will
;be sent through this piece, followed by PIPE2.PRO

;CHANGE MADE ON 01-17-09
;AS THE NEW PIPELINE WORKS FROM THE UNCOLLAPSED DATA THIS NEW VERSION
;OF PIPE1 USES THE PEFSM.FITS FILES AND CORRESPONDING PEFW.FITS FILES
;TO DO A QUICK COLLAPSE OF THE DATA BEFORE BEING FED INTO
;RADPLOT_FUNC.PRO. 

;The quick collapse is of the form
;         collapsed = sum_over_i(pefsm_i * pefw_i)/sum_over_i(pefw_i)

;Another alteration is eliminating the .rad files. Now, all the
;values are kept in the .R file with the following format.

; fiber #   Radius   dRA    dDec   Flux (normalized)   Flux (fiber median)

; The dRA and dDec values (deltaRA, deltaDec) is the offset from the
; center of the galaxy as entered by you in the procedure call.
;****************************************************************************
PRO pipe1, pointing, galRA, galDec, Data=datalist, SKIP=skipreject 
;****************************************************************************

; This code calls radplotf.pro and wgtcollapsef.pro.

;EX: 'n1600cent_D1','hh:mm:ss','+dd:mm:ss','n1600cent.list'

;POINTING: The name used in finder_code_gamma, WITHOUT the "_coords.txt"
;GALRA AND GALDEC: The galaxy RA and Dec (feed in as strings)
;DATALIST: A list that contains all the frames FOR A GIVEN POINTING. In the
;case of multiple sky-subtractions, only one needs to be sent through
;as the name.R file will hold for all of them.
;(NOTE: If the datalist keyword isn't used, the routine looks for a
;list named pointing.list.)

;NOTE THAT THEY SHOULD ALL BE THE SAME LEVEL OF SUBTRACTION AS THE
;UN-NORMALIZED FLUX IS USED AS A WEIGHT TO DETERMINE THE BINNING
;LOCATIONS IN THE VISPLOT AND OTHER ROUTINES...

;The format of the list is the same as used in pipe2.pro, so some of
;the columns are unnecessary.
;WARNING: YOU WILL NEED TO PAY ATTENTION TO THE LEVEL OF SKY-SUBTRACTION IN
;THE LIST. IF THESE LEVELS VARY MUCH OTHER ROUTINES THAT USE THE
;ESTIMATED FLUX WILL BE THROWN OFF.

;The format of the list should be as follows:
;   jm0302cc_   a   _jan08_n2.dat    jm0302.R  jm0301_  jm0303_
; The 4th column (dotR) is just the name you want to give to the
; output files.

;MODIFIED ON JULY 6, 2011: Added a switch called "SKIP". If set to
;'yes' then the fiber rejection step is turned off and the radplotF
;routine just doesn't ask for fibers to be rejected. If not used, the
;rejection step will be used.

;MODIFIED ON FEB 13, 2012: This code was slightly modified to
;accomodate the changes made to pipe2.pro which requires slightly
;different lists.

; M87: 12:30:49.4 +12:23:28
; N307: 00:56:32.6 -01:46:19
; N1600: 
;****************************************************************************
;----------------------------------------------------------------------------
nprofile = 5 ;this is the number of elements in your fiber profile

pscale = 0.005 ;decrease this number to stretch the level of the coutour plots.
             ;set > 1.0 for galaxies that are larger on the sky.
csize = 1.2 ;this sets the chararcter size in the contour plots. for the 2.7m
            ;FOV then 1.0 is fine. the HET requires a smaller one.
filesuffix = 'pefsm.fits'
computer = 'grad78'
;----------------------------------------------------------------------------
;****************************************************************************

if (n_elements(skipreject) eq 0) then skipreject = 'no'

if (n_elements(datalist) eq 0) then datalist = pointing+'.list'
readcol,datalist,format='a,x,a,a,x,x',filenames,dat,dotR
print,'The files in the list are:'
print,transpose(filenames)

file = pointing+'_coords.txt'
RA = strn(galRA) & Dec = strn(galDec)
RA = strsplit(RA,':',/EXTRACT)
Dec = strsplit(Dec,':',/EXTRACT)
RA = double((RA[0]+RA[1]/60.+RA[2]/3600.)*15)
;The + or - sign of the declination is accounted for
IF (Dec[0] LT 0) THEN Dec = double(Dec[0]-Dec[1]/60.-Dec[2]/3600.)
IF (Dec[0] GE 0) THEN Dec = double(Dec[0]+Dec[1]/60.+Dec[2]/3600.)

;The RA and Dec for each fiber is read in.
readcol,file,format='I,A,A',fib,fra,fdec
n1 = n_elements(fib)
rad = dblarr(4,n1) ;The radius array (fiber#, rad, deltaRA, deltaDec)
cor = dblarr(3,n1) ;The coordinate array (fiber#, RA, Dec)

FOR j=0,n1-1 DO BEGIN ;A loop through each fiber
    cor[0,j] = fib[j] ;dead fibers are not in this list.
    fc = fra[j]
    fd = fdec[j]
    fc = strsplit(fc,':',/EXTRACT)
    fd = strsplit(fd,':',/EXTRACT)
;The coordinates are converted to decimal degrees
    cor[1,j] = double((fc[0]+fc[1]/60.+fc[2]/3600.)*15.)
    IF (fd[0] LT 0) THEN cor[2,j] = double(fd[0]-fd[1]/60.-fd[2]/3600.)
    IF (fd[0] GT 0) THEN cor[2,j] = double(fd[0]+fd[1]/60.+fd[2]/3600.)
ENDFOR

;The radial distance from the galaxy center, deltaRA and deltaDec are
;calculated, then converted back to arcseconds.
FOR j=0,n1-1 DO BEGIN
    dx = double((cor[1,j]-RA))
    dy = double((cor[2,j]-Dec))
    rad[0,j] = uint(cor[0,j])
    rad[1,j] = sqrt(dx^2+dy^2)*3600.
    rad[2,j] = dx*3600.
    rad[3,j] = dy*3600.
;    Print, 'The radius for fiber '+strn(j+1)+ ' is: ',strn(rad[1,j])
ENDFOR

xmin = min(cor[1,*])-0.001
xmax = max(cor[1,*])+0.003
ymin = min(cor[2,*])-0.002
ymax = max(cor[2,*])+0.002

;A visual confirmation of the positions of each fiber. This should
;always look just like a plot of the IFU head.

window,0,retain=2,xsize=650,ysize=650,xpos=0,ypos=700
device,decomposed=0
loadct,0,/silent
plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
  title='Radial positions of fibers (arcsec) relative to center',$
  charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1,$
  /isotropic
loadct,4,/silent
for j=0,n1-1 do xyouts,cor[1,j],cor[2,j],strn(round(rad[1,j])),$
  color=150,charsize=csize

print,'Next ENTER deletes this plot...'
pause
wdelete,0

openw,6,pointing+'.coord'
for j=0,n1-1 do begin
    printf,6,uint(rad[0,j]),rad[1,j],rad[2,j],rad[3,j]
endfor
free_Lun,6
print,'The file '+pointing+'.coord has been written.'

;Now that the "name.coord" file is created, it is fed into a modified
;version of radplot. Here, individual fibers are rejected and the
;"name.R" files are generated.

if (skipreject eq 'no') then begin
    ans1=''
    print,'Plot contour lines? (y/n)'
    read,ans1
endif

for j=0,n_elements(filenames)-1 do begin ;a loop through each file
    data = readfits(filenames[j]+filesuffix,/silent)
    weight = readfits(filenames[j]+'pefw.fits',/silent)
    maskname = 'mask'+dat[j]

    datac = wgtcollapsef(data,weight,maskname,nprofile)
;    datac = MEDcollapseF(data,maskname)
    if (skipreject eq 'yes') then goto,skipjump
;the following is all just for plotting purposes.
    values = dblarr(n1)
    for k=0,n1-1 do values[k] = median(datac[*,k])
    nbad = where(values eq -666 or values lt 0.0)
    if (nbad[0] ne -1) then values[nbad] = 0.0
    
    out = 'y'
    repeat begin
        ival = where(values lt pscale*max(values),valct)
        if (valct gt 1) then begin
            if (max(values[ival]) ne 0) then begin
                pvalues = double(values/max(values[ival]))
            endif else begin
                print,'Out of range! That number is too small/big! Try again...'
                goto,jump1
            endelse 
        endif
        pvalues = round(pvalues*250)
        iwhite = where(pvalues ge pscale * max(pvalues),wcount)
        lvl = findgen(20)*15+min(pvalues[where(pvalues gt 0)])

        if (computer eq 'grad78') then $
        window,1,retain=2,xsize=600,ysize=550,xpos=90,ypos=700 ;for grad78
        if (computer eq 'danu') then $
        window,1,retain=2,xsize=500,ysize=450,xpos=260,ypos=700 ;for danu

        device,decomposed=0
        loadct,0,/silent
        plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
          title='Fiber Number (colored by flux level)',/isotropic,$
          charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1
        loadct,27,/silent
        if (ans1 eq 'y') then contour,pvalues,cor[1,*],cor[2,*],$
          /irregular,/overplot,levels=lvl,min_value=1,charsize=csize
        for k=0,n1-1 do xyouts,cor[1,k],cor[2,k],strn(uint(cor[0,k])),color=pvalues[k],$
          charsize=csize
        if (wcount ne 0) then begin
            loadct,0,/silent
            for k=0,wcount-1 do xyouts,cor[1,iwhite[k]],$
              cor[2,iwhite[k]],strn(iwhite[k]+1),charsize=csize
        endif

        if (computer eq 'grad78') then $
        window,2,retain=2,xsize=600,ysize=550,xpos=690,ypos=700 ;for grad78
        if (computer eq 'danu') then $
        window,2,retain=2,xsize=500,ysize=450,xpos=770,ypos=700 ;for danu

        device,decomposed=0
        loadct,0,/silent
        plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
          title='Relative Flux in Each Fiber for '+filenames[j],/isotropic,$
          charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1
        loadct,27,/silent
        if (ans1 eq 'y') then contour,pvalues,cor[1,*],cor[2,*],$
          /irregular,/overplot,levels=lvl,min_value=1,charsize=csize
        for k=0,n1-1 do xyouts,cor[1,k],cor[2,k],strn(pvalues[k]),$
          color=pvalues[k],charsize=csize
        if (wcount ne 0) then begin
            loadct,0,/silent
            for k=0,wcount-1 do xyouts,cor[1,iwhite[k]],$
              cor[2,iwhite[k]],strn(pvalues[iwhite[k]]),charsize=csize
        endif
        jump1:
        if (skipreject eq 'no') then begin
            print,'Adjust the contrast? (y/n):'
            read,out
            if (out eq 'y') then begin
                print,'The current pscale value is '+strn(pscale)
                print,'Lower this value to stretch the constrast...'
                print,'Raise this value to shrink the contrast...'
                print,'Enter a new pscale value:'
                read,pscale
            endif
        endif else out = 'n'
    endrep until (out eq 'n')

; ps of the fiber positions is made...
    set_plot,'ps'
    device,filename=dotR[j]+'_fibermap.ps',/color
    loadct,0,/silent
    plot,cor[1,*],cor[2,*],psym=3,xtitle='RA',ytitle='Dec',/nodata,$
      title='VP2 Fiber Map for '+dotR[j],xthick=3,ythick=3,charthick=3,$
      charsize=1.2,xrange=[xmin,xmax],yrange=[ymin,ymax],xstyle=1,ystyle=1
    loadct,27,/silent
    for k=0,n1-1 do xyouts,cor[1,k],cor[2,k],strn(k+1),charthick=3,$
      color=pvalues[k]+1,charsize=csize-0.2
    device,/close_file
    set_plot,'x'

    skipjump:
;Radplot is now a function. The string can be either 'm' for median
;or 't' for total. This is simply what radplot.pro plots to aid in the
;visual rejection. i.e. the median or total of the fiber vs. radius 
    nameout = radplotf(pointing+'.coord',datac,dotR[j],'m',skip=skipreject)

endfor

print,''
print,'PIPE1 has successfully finished'
print,''

wait,3.0
while !d.window ne -1 do wdelete, !d.window

stop
END

