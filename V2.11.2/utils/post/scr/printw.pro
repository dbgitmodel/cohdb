pro printw, visfile, device=wdevice, outname=outname, printer=wprinter, $
            date=date, nobatch=nobatch, landscape=wland, $
            color=colortable
;**********************************************************************
;
;    *printw* 
;
;       Author - patrick luyten
;
;       Last update - 3 Jun 2003       @(#)printw.pro
;
;       Description - IDL procedure to plot a succession of
;                     muplot/mumap2 files
;
;**********************************************************************
;
;*    Arguments
;
;----------------------------------------------------------------------
;     Name      Type      Purpose
;     ----      ----      -------
;    *visfile*  string    file with names of parameter files
;    *device*   string    type of output device
;    *outname*  string    name of parameter output file
;    *printer*  string    printer type
;    *date*     integer   switch to select date on the image (mumap2 only)
;    *nobatch*  integer   switch to disable/enable batch mode
;    *landscape*integer  switch to select landscape (on) or portrait (off)
;                         (mumap2 only)
;    *color*    int or string specifies color table or file
;
;---------------------------------------------------------------------------

; initialise parameters
visline =string(0)
viscase=string(0)
plotprog = string(0)
toplot = string(0)
plotnam = string(0)
wparfile = string(0)

; treatment of input-line parameters
if (n_elements(visfile) eq 0) then visfile = 'files.vis'
if (n_elements(wdevice) eq 0) then wdevice = 'screen'
if (n_elements(outname) eq 0) then $
  case wdevice of
    'pict' :    outname = 'mu-image.pict'
    'tv'   :    outname = 'mu-image.tv'
    'tv8'  :    outname = 'mu-image.tv8'
    't24'  :    outname = 'mu-image.tv24'
    'ppm'  :    outname = 'mu-image.ppm'
    'pscolor' : outname = 'POSTS.COLOR'
    'eps'  :    outname = 'EPS'
    'epscolor': outname = 'EPS.COLOR'
    'gif':      ourtname = 'gif'
    else :      outname = 'POSTS'
  endcase
if (n_elements(wprinter) eq 0) then wprinter = 'screen'
if (n_elements(date) eq 0) then nodate = 1
if (n_elements(nobatch) eq 0) then batch = 1
if (n_elements(wland) eq 0) then wland = 0
if (n_elements(colortable) eq 0) then colortable = -1

; check arguments
if  ((wdevice ne 'screen') and (wdevice ne 'ps') and (wdevice ne 'pscolor') $
 and (wdevice ne 'eps') and (wdevice ne 'epscolor') $
 and (wdevice ne 'tek4107') and (wdevice ne 'pict') and (wdevice ne 'tv') $
 and (wdevice ne 'tv8') and (wdevice ne 'tv24') and (wdevice ne 'ppm') $
 and (wdevice ne 'gif')) $
then begin
  print, 'Wrong output device : ', wdevice
  stop
endif
if  ((wprinter ne 'screen') and (wprinter ne 'RICOH_B') and (wprinter ne 'RICOH_C')) then begin
  print, 'Wrong printer type : ', wprinter
  stop
endif
if  ((wdevice ne 'ps') and (wdevice ne 'pscolor') and (wdevice ne 'eps') $
 and (wdevice ne 'epscolor')) then wprinter = 'screen'

if (wprinter eq 'screen') then begin
  print, 'Plotting on screen'
endif else begin
  print, 'Printing on ', wprinter
endelse

; open 'files.vis'
openr, 99, visfile
next_test :
repeat begin
  readf, 99, viscase  
  viscase = strcompress(viscase)
  if (strupcase(viscase) eq 'END') then goto, end_vis
  ipos = strpos(viscase,'!')
endrep until (ipos ne 0)

; open 'param.vis'
openr, 98, viscase

; read plot information
while (not (eof(98) or (strupcase(plotprog) eq 'END'))) do begin
   readf, 98, visline
   visline = strtrim(visline,2)
   ispace = strpos(visline,' ')
   if (ispace gt 0) then plotprog=strmid(visline,0,ispace) else goto, quit
   visline = strtrim(strmid(visline,ispace,strlen(visline)),2)
   ispace = strpos(visline,' ')
   wparfile = strmid(visline,0,ispace)
   visline = strtrim(strmid(visline,ispace,strlen(visline)),2)
   ispace = strpos(visline,' ')
   toplot = strmid(visline,0,ispace)
   plotnam = strtrim(strmid(visline,ispace,strlen(visline)),2)
   print, plotprog,'-',wparfile,'-',toplot,'-',plotnam

; print figures
   if (toplot eq 'y') then begin
     if (wdevice eq 'screen') then begin
       print,'hit q for quit, n to skip this plot, any other key to show plot'
       anykey = get_kbrd(1)
       if (strlowcase(anykey) eq 'q') then goto, quit
       if (strlowcase(anykey) eq 'n') then goto, next_plot
     endif
     if (plotprog eq 'muplot') then $
       muplot, wparfile, landscape=wland, device=wdevice, $
                 outname=outname, printer=wprinter
     if (plotprog eq 'mumap2') then begin
       ctable = colortable
       mumap2, wparfile, landscape=wland, nodate = nodate, $
               batch=batch, device=wdevice, outname=outname, color=ctable
       if (wprinter ne 'screen') then begin
         cmd = 'lpr -P'+wprinter+' '+outname
         spawn,cmd
       endif
      endif
   endif
next_plot :
endwhile
quit:
  close, 98
if (not(eof(99))) then goto, next_test

end_vis :
  close, 99

end
