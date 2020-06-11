pro muplot, filename, landscape=wland,device=wdevice, $
              outname=outname, printer=wprinter
;**********************************************************************
;
;    *muplotplot*      IDL PLOTTING PROCEDURE
;
;	AUTHOR - PATRICK LUYTEN
;
;       LAST UPDATE - 3 Jun 2003       @(#)muplot.pro
;
;	DESCRIPTION - READS PARAMETERS FROM INPUT FILE filename
;		
;**********************************************************************
;
;*    ARGUMENTS
;
;----------------------------------------------------------------------
;     NAME      TYPE      PURPOSE
;     ----      ----      -------
;    *filename* string    name of parameter file
;    *device*   string    type of output device
;    *outname*  string    name of output file
;    *printer*  string    printer type
;
;----------------------------------------------------------------------

; initialise parameters
minval = 0.0
nrecmax = 10000
ncurves = 0
datfil0 = string(0)
xtit = string(0)
ytit = string(0)
maintit = string(0)
subtit = string(0)
ixstyle = 0
ixrange = 0
xdatmin = 0.0
xdatmax = 0.0
iystyle = 0
iyrange = 0
ydatmin = 0.0
ydatmax = 0.0

; treatment of input parameters
if (n_elements(filename) eq 0) then filename='muplot.param'
if (n_elements(wdevice) eq 0) then wdevice='screen'
if (n_elements(wland) eq 0) then wland = 0
if  ((wdevice ne 'ps') $
 and (wdevice ne 'pscolor') and (wdevice ne 'eps') $
 and (wdevice ne 'epscolor')) then wdevice='screen'

;if (landscape eq 1) then wland = 1 else wland = 0

if (n_elements(outname) eq 0) then $
  case wdevice of
    'pscolor' : outname = 'POSTS.COLOR'
    'eps'     : outname = 'EPS'
    'epscolor': outname = 'EPS.COLOR'
    else      : outname = 'POSTS'
  endcase
if (n_elements(wprinter) eq 0) then wprinter='screen'

;check arguments
if  ((wprinter ne 'screen') and (wprinter ne 'ricoh4502')) then begin
  print, 'Wrong printer type : ', wprinter
  stop
endif

; open parameter file
iunit = 97
openr, iunit, filename

; read number of curves
readf, iunit, ncurves, minval
if (minval gt 0.0) then begin
  minval = -minval
endif
minval = 0.99*minval

; initialise arrays
datfiles = strarr(ncurves)
xdat = fltarr(ncurves,nrecmax)
ydat = fltarr(ncurves,nrecmax)
xmin = fltarr(ncurves)
ymin = fltarr(ncurves)
xmax = fltarr(ncurves)
ymax = fltarr(ncurves)
linepsyms = intarr(ncurves)
linest = intarr(ncurves)
psymst = intarr(ncurves)
nrec = intarr(ncurves)

; read data file names
for n=0,ncurves-1 do begin
   readf, iunit, datfil0
   print, '  data file: ', datfil0
   datfiles(n) = datfil0
endfor

; read titles
readf, iunit, xtit 
readf, iunit, ytit
readf, iunit, maintit 
readf, iunit, subtit

; read xrange and yrange parameters
readf, iunit, ixstyle, ixrange, xdatmin, xdatmax
readf, iunit, iystyle, iyrange, ydatmin, ydatmax

; read linestyles
readf, iunit, linepsyms
for n=0,ncurves-1 do begin
  psymst(n) = linepsyms(n)/10
  linest(n) = linepsyms(n) - 10*psymst(n) - 1
endfor

; close parameter file
close, iunit

; read data
for n=0,ncurves-1 do begin
  openr, iunit, strtrim(datfiles(n))
  irec = 0   
  while not eof(iunit) do begin
    readf, iunit, xdat0, ydat0
    xdat(n,irec) = xdat0 & ydat(n,irec) = ydat0
    if (irec lt nrecmax) then irec = irec+1 else stop, 'too many records'
  endwhile
  nrec(n) = irec
  close, iunit
endfor

; titles of the plot
maintit = strtrim(maintit,2)
subtit = strtrim(subtit,2)
xtit = strtrim(xtit,2)
ytit = strtrim(ytit,2)

; xrange and yrange parameters
if (ixrange eq 0) then begin
  xdatmin = 0.0
  xdatmax = 0.0
endif else if (ixrange eq 1) then begin
  for n=0,ncurves-1 do begin
    xmin(n) = min(xdat(n,0:nrec(n)-1))
    xmax(n) = max(xdat(n,0:nrec(n)-1))
  endfor
  xdatmin = min(xmin)
  xdatmax = max(xmax)
endif
if (iyrange eq 0) then begin
  ydatmin = 0.0
  ydatmax = 0.0
endif else if (iyrange eq 1) then begin
  for n=0,ncurves-1 do begin
    ymin(n) = min(ydat(n,0:nrec(n)-1))
    ymax(n) = max(ydat(n,0:nrec(n)-1))
  endfor
  ydatmin = min(ymin)
  ydatmax = max(ymax)
endif

; set device 
olddevice = !D.name
if (wdevice eq 'screen') then begin
  set_plot, 'X'
endif else begin
  set_plot, 'ps'
endelse
case wdevice of
  'srcreen' :  device, landscape=wland
  'ps' :       device, encapsulated=0, landscape=wland, filename=outname
  'pscolor' :  device, encapsulated=0, landscape=wland, /color, filename=outname
  'eps' :      device, /encapsulated,  landscape=wland, filename=outname
  'epscolor' : device, /encapsulated,  landscape=wland, /color, filename=outname
  else :
endcase

if (linest(0) gt -1) then begin
  plot, xdat(0,0:nrec(0)-1), ydat(0,0:(nrec(0)-1)), xstyle=ixstyle, $
    ystyle=iystyle, xrange=[xdatmin,xdatmax], yrange=[ydatmin,ydatmax], $
    min_value=minval,linestyle=linest(0), title=maintit, subtitle=subtit, $
    xtitle=xtit, ytitle=ytit
  if (psymst(0) gt 0) then begin
     oplot, xdat(0,0:nrec(0)-1), ydat(0,0:(nrec(0)-1)), min_value=minval, $
       psym = psymst(0)
  endif
endif else begin
  plot, xdat(0,0:nrec(0)-1), ydat(0,0:(nrec(0)-1)), xstyle=ixstyle, $
    ystyle=iystyle, xrange=[xdatmin,xdatmax], yrange=[ydatmin,ydatmax], $
    min_value=minval,psym=psymst(0), title=maintit, subtitle=subtit, $
    xtitle=xtit, ytitle=ytit
endelse

for n=1,ncurves-1 do begin
  if (linest(n) gt -1) then begin
    oplot, xdat(n,0:nrec(n)-1), ydat(n,0:nrec(n)-1), min_value=minval, $
      linestyle=linest(n)
    if (psymst(n) gt 0) then begin
      oplot, xdat(n,0:nrec(n)-1), ydat(n,0:nrec(n)-1), min_value=minval, $
        psym=psymst(n)
    endif
  endif else begin
      oplot, xdat(n,0:nrec(n)-1), ydat(n,0:nrec(n)-1), min_value=minval, $
        psym=psymst(n)
  endelse
endfor

; close device and send to printer (if necessary)
if (wdevice ne 'screen') then device,/close
if (wprinter ne 'screen') then begin
   cmd = 'lpr -P'+wprinter+' '+outname
   spawn, cmd
endif

set_plot, olddevice

return

end
