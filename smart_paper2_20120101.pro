pro smart_paper2_20120101, res1tore=res1tore, res2tore=res2tore, res3store=res3store, res4store=res4store, $
	skiptoplot2=skiptoplot2

datapath='~/science/papers/active_regions_2_cycle/data/'
plotpath='~/science/papers/active_regions_2_cycle/images/20110330/'

;STRUCTURE OF ALL SMART DETECTIONS-------------------------------------------->
if not keyword_set(res1tore) then begin
	;Generate structure of all AR detections
	ffdetect1='~/science/bradford/sipwork_paper/flare_associate/flare_assoc_structures_1996-2002_2009-2011.sav'
	ffdetect2='~/science/bradford/sipwork_paper/flare_associate/flare_assoc_structures_2003-2008.sav'
	restore,ffdetect1,/ver
	
	;combine 2 data sets
	arstr1=ARSTR_ARR
	yrs=strmid(arstr_arr.time,7,4)
	w2002end=max(where(yrs eq '2002'))
	
	restore,ffdetect2,/ver
	ARSTR_ARR=[arstr1[0l:w2002end], ARSTR_ARR, arstr1[w2002end+1l:*]]
	
	;fine where region boundaries between +- 60 degrees lon
	wbtwn60=where((ARSTR_ARR.extstr.hglon)[0,*] ge -60. and (ARSTR_ARR.extstr.hglon)[1,*] le 60.)
	
	;create time array
	ARSTR_ARR_tim=anytim(ARSTR_ARR.time)
	
	save,ARSTR_ARR,ARSTR_ARR_tim,wbtwn60,file=datapath+'smart_detect_all_20110330.sav'
endif else restore, datapath+'smart_detect_all_20110330.sav',/ver

;---take only detections within 60 degrees---
ARSTR_ARR=ARSTR_ARR[wbtwn60]
ARSTR_ARR_tim=ARSTR_ARR_tim[wbtwn60]

;setup plot colors and symbols
setcolors,/sys,/silent,/quiet
pcolor=!black
ncolor=!red
scolor=!blue
nsym=-1
ssym=-4

;CONSTANTS
secperyr=3600.*24.*365.

;CHECK FOR PLOT 2 SKIP KEY WORD
if keyword_set(skiptoplot2) then begin & print,'SKIPPING TO PLOT 2...........' & goto,gotoplot2 & endif

;SMART TIME SERIES (FLUX, FLUX NORTH, FLUX SOUTH, NET FLUX NORTH SOUTH)
if not keyword_set(res2tore) then begin
	restore,'~/science/data/cycle23_sav/sunspot_monthly_num_sdic.sav'; sstim,ssnum
	tseriesbin=3600.*27.*24. ;bin by 27 days
	;FLUX
	smart_bin_time, arstr_arr.bflux, ARSTR_ARR_tim, flux27dseries, tflux27dseries, bin=tseriesbin,/total
	wnorthdetect=where(arstr_arr.hglat ge 0)
	wsouthdetect=where(arstr_arr.hglat lt 0)
	smart_bin_time, arstr_arr[wnorthdetect].bflux, ARSTR_ARR_tim[wnorthdetect], fluxN27dseries, tfluxN27dseries, bin=tseriesbin,/total
	smart_bin_time, arstr_arr[wsouthdetect].bflux, ARSTR_ARR_tim[wsouthdetect], fluxS27dseries, tfluxS27dseries, bin=tseriesbin,/total
	;NET FLUX
	smart_bin_time, arstr_arr[wnorthdetect].bfluxpos-abs(arstr_arr[wnorthdetect].bfluxneg), ARSTR_ARR_tim[wnorthdetect], netfluxN27dseries, tnetfluxN27dseries, bin=tseriesbin,/total
	smart_bin_time, arstr_arr[wsouthdetect].bfluxpos-abs(arstr_arr[wsouthdetect].bfluxneg), ARSTR_ARR_tim[wsouthdetect], netfluxS27dseries, tnetfluxS27dseries, bin=tseriesbin,/total
	save,sstim,ssnum,tseriesbin,flux27dseries,tflux27dseries, $
		fluxN27dseries,tfluxN27dseries,fluxS27dseries,tfluxS27dseries,netfluxN27dseries, $
		tnetfluxN27dseries, netfluxS27dseries,tnetfluxS27dseries, $
		file=datapath+'smart_timeseries_20110330.sav'
endif else restore,datapath+'smart_timeseries_20110330.sav',/ver

;SMART BUTTERFLY DIAGRAMS----------------------------------------------------->
if not keyword_set(res3tore) then begin
	;bin by 27 days and 1 degree latitude
	ybin=1. & xbin=27.*3600.*24.
	xran=anytim(['1-jun-1996','1-jan-2011'])
	yran=[-60,60]
	
	;Generate butterfly diagrams from SMART detections (total flux and flux imbalance)
	xx=smart_bin2d_simple(xbin=xbin,ybin=ybin,xran=xran,yran=yran,struct=arstr_arr, $ ;,/years
		/latitude,timarr=ARSTR_ARR_tim,path=datapath,outfile=outfile,filemod='_27day_bin')
	restore,outfile,/ver
endif else restore,'~/science/papers/active_regions_2_cycle/data/butterfly_diags_lat_27day_bin.sav'


;---PLOT 0. - solar cycle time series versus flux----------------------------->
setplotenv,xs=13,ys=13,/ps,file=plotpath+'plot_0_cycletseries_vs_fluxbflydiag.eps'
setcolors,/sys,/silent,/quiet
	mintimseries=min(tflux27dseries)
	postseries=[.13,.525,.95,.95]
	tseriesrange=[min(xarr),max(xarr)]-mintimseries;[anytim('1-jun-1996'),anytim('1-jan-2011')]-mintimseries
	;plot time series here: north total flux on disk per day, south total flux on disk per day
	utplot,tflux27dseries-mintimseries,flux27dseries*1d16,mintimseries,ps=10,color=pcolor,xtit='',xtickname=strarr(10)+' ',xran=tseriesrange,/xsty,pos=postseries,ytit=textoidl('\Phi_{TOT}')+' [Mx]'
	oplot,sstim-mintimseries,ssnum/max(ssnum)*max(flux27dseries)*1d16,color=!gray
	oplot,tfluxN27dseries-mintimseries,fluxN27dseries*1d16,ps=nsym,color=ncolor
	oplot,tfluxS27dseries-mintimseries,fluxS27dseries*1d16,ps=ssym,color=scolor
	
	legend,['Total', 'North Hemi.', 'South Hemi.', 'Sunspot #'],/top,/right,color=[0,ncolor,scolor,!gray],psym=[0,nsym,ssym,0]
loadct,0
	;flux butterfly diag
	posfluximg=[.13,.1,.95,.525]
	plot_image,(-1.)*FLUXIMAGE^(.4) > (-300.),pos=posfluximg,color=255,/noerase
	utplot,xarr-mintimseries,yarr,mintimseries,/nodat,pos=posfluximg,/noerase,/xsty,/ysty,ytit='Latitude',xtit='Year (Begins 1 June 1996)'
	
	cbchars=2
;	colorbar,maxrange=(300.)^(1./.4)*1d16,minrange=0.,/invert,pos=[.6,.14,.8,.16],color=0,chars=!p.charsize
	colimg=(-1.)*(findgen(256)/255.*(1.56d6-min(FLUXIMAGE))+min(FLUXIMAGE))^(.4) > (-300.)
	plot_image,transpose(rebin(reform(colimg,1,256),10,256,/samp)),pos=[.67,.47,.91,.49],/noerase,xticklen=.0001,yticklen=.0001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' ';,/nosq,xsty=4,ysty=4
	plot, findgen(10),/nodat,/noerase,/xsty,xran=[0,1],/ysty,pos=[.67,.47,.91,.49], $
		xticklen=.0001,yticklen=.0001,xtit='Flux [Mx]',xtickname=strarr(10)+' ',ytickname=strarr(10)+' ',chars=cbchars;,/nosq,xticks=3,
	vline,[.333,.837295],col=255;,(1.-.837295)]
	xyouts,0,-1.5,'0        1E21           1E22',/data,chars=cbchars

closeplotenv
;---END PLOT 0.--------------------------------------------------------------->
stop

;ANALYZE BUTTERFLY DIAGRAMS--------------------------------------------------->
smthrad=[0,0];[3,3]
missing=-10000.
nFLUXIMAGE=FLUXIMAGE[*,60:*];smooth(FLUXIMAGE[*,60:*],smthrad)
nFLUXgt22IMAGE=FLUXgt22IMAGE[*,60:*];smooth(FLUXgt22IMAGE[*,60:*],smthrad)
nFLUXgt5d22IMAGE=FLUXgt5d22IMAGE[*,60:*];smooth(FLUXgt5d22IMAGE[*,60:*],smthrad)
nFLUXgt23IMAGE=FLUXgt23IMAGE[*,60:*];smooth(FLUXgt23IMAGE[*,60:*],smthrad)
sFLUXIMAGE=FLUXIMAGE[*,0:59];smooth(FLUXIMAGE[*,0:59],smthrad)
sFLUXgt22IMAGE=FLUXgt22IMAGE[*,0:59];smooth(FLUXgt22IMAGE[*,0:59],smthrad)
sFLUXgt5d22IMAGE=FLUXgt5d22IMAGE[*,0:59];smooth(FLUXgt5d22IMAGE[*,0:59],smthrad)
sFLUXgt23IMAGE=FLUXgt23IMAGE[*,0:59];smooth(FLUXgt23IMAGE[*,0:59],smthrad)

;Fill missing values in maps
wmiss=where(total(FLUXIMAGE,2) le 0) & wmiss=wmiss[where(wmiss lt 150)] ;use nFLUXIMAGE?
nFLUXIMAGE[wmiss,*]=missing
nFLUXgt22IMAGE[wmiss,*]=missing
nFLUXgt5d22IMAGE[wmiss,*]=missing
nFLUXgt23IMAGE[wmiss,*]=missing
sFLUXIMAGE[wmiss,*]=missing
sFLUXgt22IMAGE[wmiss,*]=missing
sFLUXgt5d22IMAGE[wmiss,*]=missing
sFLUXgt23IMAGE[wmiss,*]=missing
FILL_MISSING, nFLUXIMAGE, missing, 1
FILL_MISSING, nFLUXgt22IMAGE, missing, 1
FILL_MISSING, nFLUXgt5d22IMAGE, missing, 1
FILL_MISSING, nFLUXgt23IMAGE, missing, 1
FILL_MISSING, sFLUXIMAGE, missing, 1
FILL_MISSING, sFLUXgt22IMAGE, missing, 1
FILL_MISSING, sFLUXgt5d22IMAGE, missing, 1
FILL_MISSING, sFLUXgt23IMAGE, missing, 1

;make contour of cycle 23
fluxlevel=0
nFLUXcont=fltarr((size(nFLUXIMAGE))[1],(size(nFLUXIMAGE))[2])
nFLUXgt22cont=nFLUXcont
nFLUXgt5d22cont=nFLUXcont
nFLUXgt23cont=nFLUXcont
sFLUXcont=nFLUXcont
sFLUXgt22cont=nFLUXcont
sFLUXgt5d22cont=nFLUXcont
sFLUXgt23cont=nFLUXcont
nFLUXcont[where(nFLUXIMAGE gt fluxlevel)]=1.
nFLUXgt22cont[where(nFLUXgt22IMAGE gt fluxlevel)]=1.
nFLUXgt5d22cont[where(nFLUXgt5d22IMAGE gt fluxlevel)]=1.
nFLUXgt23cont[where(nFLUXgt23IMAGE gt fluxlevel)]=1.
sFLUXcont[where(sFLUXIMAGE gt fluxlevel)]=1.
sFLUXgt22cont[where(sFLUXgt22IMAGE gt fluxlevel)]=1.
sFLUXgt5d22cont[where(sFLUXgt5d22IMAGE gt fluxlevel)]=1.
sFLUXgt23cont[where(sFLUXgt23IMAGE gt fluxlevel)]=1.

;Find largest contour
nFLUXcontsep=smart_cont_sep(nFLUXcont, contlevel=.5, vthresh=1, areathresh=1)
nFLUXcontsep[where(nFLUXcontsep ne 1)]=0
sFLUXcontsep=smart_cont_sep(sFLUXcont, contlevel=.5, vthresh=1, areathresh=1)
sFLUXcontsep[where(sFLUXcontsep ne 1)]=0

;find centroid locations
xyrcoord,size(nFLUXIMAGE),dum1,nyy
nfluxcent=total(nFLUXIMAGE*nFLUXcontsep*nyy,2)/total(nFLUXIMAGE*nFLUXcontsep,2)
nfluxgt22cent=total(nFLUXgt22IMAGE*nFLUXcontsep*nyy,2)/total(nFLUXgt22IMAGE*nFLUXcontsep,2)
nfluxgt5d22cent=total(nFLUXgt5d22IMAGE*nFLUXcontsep*nyy,2)/total(nFLUXgt5d22IMAGE*nFLUXcontsep,2)
nfluxgt23cent=total(nFLUXgt23IMAGE*nFLUXcontsep*nyy,2)/total(nFLUXgt23IMAGE*nFLUXcontsep,2)
sfluxcent=total(sFLUXIMAGE*sFLUXcontsep*rot(-nyy,180),2)/total(sFLUXIMAGE*sFLUXcontsep,2)
sfluxgt22cent=total(sFLUXgt22IMAGE*sFLUXcontsep*rot(-nyy,180),2)/total(sFLUXgt22IMAGE*sFLUXcontsep,2)
sfluxgt5d22cent=total(sFLUXgt5d22IMAGE*sFLUXcontsep*rot(-nyy,180),2)/total(sFLUXgt5d22IMAGE*sFLUXcontsep,2)
sfluxgt23cent=total(sFLUXgt23IMAGE*sFLUXcontsep*rot(-nyy,180),2)/total(sFLUXgt23IMAGE*sFLUXcontsep,2)

;find min max edge of each map
nfluxmm=fltarr(2,(size(nFLUXIMAGE))[1])-10000.
nfluxmmgt22=nfluxmm
nfluxmmgt5d22=nfluxmm
nfluxmmgt23=nfluxmm
sfluxmm=nfluxmm
sfluxmmgt22=nfluxmm
sfluxmmgt5d22=nfluxmm
sfluxmmgt23=nfluxmm

for i=0,(size(nFLUXIMAGE))[1]-1 do begin
	wthis=where((nFLUXcont*nFLUXcontsep)[i,*] gt 0)
	wthisgt22=where((nFLUXgt22cont*nFLUXcontsep)[i,*] gt 0)
	wthisgt5d22=where((nFLUXgt5d22cont*nFLUXcontsep)[i,*] gt 0)
	wthisgt23=where((nFLUXgt23cont*nFLUXcontsep)[i,*] gt 0)
	if wthis[0] ne -1 then nfluxmm[*,i]=minmax(wthis)
	if wthisgt22[0] ne -1 then nfluxmmgt22[*,i]=minmax(wthisgt22)
	if wthisgt5d22[0] ne -1 then nfluxmmgt5d22[*,i]=minmax(wthisgt5d22)
	if wthisgt23[0] ne -1 then nfluxmmgt23[*,i]=minmax(wthisgt23)
	swthis=where((sFLUXcont*sFLUXcontsep)[i,*] gt 0)
	swthisgt22=where((sFLUXgt22cont*sFLUXcontsep)[i,*] gt 0)
	swthisgt5d22=where((sFLUXgt5d22cont*sFLUXcontsep)[i,*] gt 0)
	swthisgt23=where((sFLUXgt23cont*sFLUXcontsep)[i,*] gt 0)
	if swthis[0] ne -1 then sfluxmm[*,i]=minmax(swthis)
	if swthisgt22[0] ne -1 then sfluxmmgt22[*,i]=minmax(swthisgt22)
	if swthisgt5d22[0] ne -1 then sfluxmmgt5d22[*,i]=minmax(swthisgt5d22)
	if swthisgt23[0] ne -1 then sfluxmmgt23[*,i]=minmax(swthisgt23)
endfor

nfluxmm[where(nfluxmm eq missing)]=0./0.
nfluxmmgt22[where(nfluxmmgt22 eq missing)]=0./0.
nfluxmmgt5d22[where(nfluxmmgt5d22 eq missing)]=0./0.
nfluxmmgt23[where(nfluxmmgt23 eq missing)]=0./0.
sfluxmm[where(sfluxmm eq missing)]=0./0.
sfluxmmgt22[where(sfluxmmgt22 eq missing)]=0./0.
sfluxmmgt5d22[where(sfluxmmgt5d22 eq missing)]=0./0.
sfluxmmgt23[where(sfluxmmgt23 eq missing)]=0./0.

;set limits for fitting
nmin=min(where(finite(nfluxcent) eq 1)) & nlim=47 & print,'N limit equator edge= '+anytim(xarr[nmin],/vms)+' '+anytim(xarr[nlim],/vms)
smin=min(where(finite(sfluxcent) eq 1)) & slim=56 & print,'S limit equator edge= '+anytim(xarr[smin],/vms)+' '+anytim(xarr[slim],/vms)
ncentlim=max(where(finite(nfluxcent[0:*]) eq 1))
scentlim=max(where(finite(sfluxcent[0:*]) eq 1))

;Find where each line to fit is finite
wncent=where(finite(nfluxcent[0:*]) eq 1)
wnlo=where(finite(nfluxmm[0,0:nlim]) eq 1)
wnhi=where(finite(nfluxmm[1,0:*]) eq 1)
wnlogt22=where(finite(nfluxmmgt22[0,0:nlim]) eq 1)
wnhigt22=where(finite(nfluxmmgt22[1,0:*]) eq 1)
wnlogt5d22=where(finite(nfluxmmgt5d22[0,0:nlim]) eq 1)
wnhigt5d22=where(finite(nfluxmmgt5d22[1,0:*]) eq 1)
wscent=where(finite(sfluxcent[0:*]) eq 1)
wslo=where(finite(sfluxmm[0,0:*]) eq 1)
wshi=where(finite(sfluxmm[1,0:slim]) eq 1)
wslogt22=where(finite(sfluxmmgt22[0,0:*]) eq 1)
wshigt22=where(finite(sfluxmmgt22[1,0:slim]) eq 1)
wslogt5d22=where(finite(sfluxmmgt5d22[0,0:*]) eq 1)
wshigt5d22=where(finite(sfluxmmgt5d22[1,0:slim]) eq 1)

;Fit lines to each contour edge
nlincent=linfit((xarr[0:*])[wncent]/secperyr,(nfluxcent[0:*])[wncent], MEASURE=fltarr(n_elements(wncent))+.5, chisq=ncentchi)
nlinlo=linfit((xarr[0:nlim])[wnlo]/secperyr,(nfluxmm[0,0:nlim])[wnlo], MEASURE=fltarr(n_elements(wnlo))+.5, chisq=nlochi)
nlinhi=linfit((xarr[0:*])[wnhi]/secperyr,(nfluxmm[1,0:*])[wnhi], MEASURE=fltarr(n_elements(wnhi))+.5, chisq=nhichi)
nlinlogt22=linfit((xarr[0:nlim])[wnlogt22]/secperyr,(nfluxmmgt22[0,0:nlim])[wnlogt22], MEASURE=fltarr(n_elements(wnlogt22))+.5, chisq=nlogt22chi)
nlinhigt22=linfit((xarr[0:*])[wnhigt22]/secperyr,(nfluxmmgt22[1,0:*])[wnhigt22], MEASURE=fltarr(n_elements(wnhigt22))+.5, chisq=nhigt22chi)
nlinlogt5d22=linfit((xarr[0:nlim])[wnlogt5d22]/secperyr,(nfluxmmgt5d22[0,0:nlim])[wnlogt5d22], MEASURE=fltarr(n_elements(wnlogt5d22))+.5, chisq=nlogt5d22chi)
nlinhigt5d22=linfit((xarr[0:*])[wnhigt5d22]/secperyr,(nfluxmmgt5d22[1,0:*])[wnhigt5d22], MEASURE=fltarr(n_elements(wnhigt5d22))+.5, chisq=nhigt5d22chi)
slincent=linfit((xarr[0:*])[wscent]/secperyr,(sfluxcent[0:*])[wscent], MEASURE=fltarr(n_elements(wscent))+.5, chisq=scentchi)
slinlo=linfit((xarr[0:slim])[wslo]/secperyr,(sfluxmm[0,0:slim])[wslo]-60., MEASURE=fltarr(n_elements(wslo))+.5, chisq=slochi)
slinhi=linfit((xarr[0:*])[wshi]/secperyr,(sfluxmm[1,0:slim])[wshi]-60., MEASURE=fltarr(n_elements(wshi))+.5, chisq=shichi)
slinlogt22=linfit((xarr[0:slim])[wslogt22]/secperyr,(sfluxmmgt22[0,0:slim])[wslogt22]-60., MEASURE=fltarr(n_elements(wslogt22))+.5, chisq=slogt22chi)
slinhigt22=linfit((xarr[0:*])[wshigt22]/secperyr,(sfluxmmgt22[1,0:slim])[wshigt22]-60., MEASURE=fltarr(n_elements(wshigt22))+.5, chisq=shigt22chi)
slinlogt5d22=linfit((xarr[0:slim])[wslogt5d22]/secperyr,(sfluxmmgt5d22[0,0:slim])[wslogt5d22]-60., MEASURE=fltarr(n_elements(wslogt5d22))+.5, chisq=slogt5d22chi)
slinhigt5d22=linfit((xarr[0:*])[wshigt5d22]/secperyr,(sfluxmmgt5d22[1,0:slim])[wshigt5d22]-60., MEASURE=fltarr(n_elements(wshigt5d22))+.5, chisq=shigt5d22chi)

;use MP to get uncertainties on parameters
pnlincent=mpfitexpr('P[0]*X+P[1]', (xarr[0:*])[wncent]/secperyr, (nfluxcent[0:*])[wncent], fltarr(n_elements(wncent))+0.5, perror=nlincenterr)
pnlinlo=mpfitexpr('P[0]*X+P[1]', (xarr[0:nlim])[wnlo]/secperyr,(nfluxmm[0,0:nlim])[wnlo], fltarr(n_elements(wnlo))+0.5, perror=nlinloerr)
pslincent=mpfitexpr('P[0]*X+P[1]', (xarr[0:*])[wscent]/secperyr,(sfluxcent[0:*])[wscent], fltarr(n_elements(wscent))+0.5, perror=slincenterr)
pslinhi=mpfitexpr('P[0]*X+P[1]', (xarr[0:*])[wshi]/secperyr,(sfluxmm[1,0:slim])[wshi]-60., fltarr(n_elements(wshi))+0.5, perror=slinhierr)
;print,pnlincent
;print,pnlinlo  
;print,pslincent
;print,pslinhi
;print,nlincenterr
;print,nlinloerr
;print,slincenterr
;print,slinhierr

save,nlincent,nlinlo,nlinhi,nlinlogt22,nlinhigt22,nlinlogt5d22,nlinhigt5d22,slincent,slinlo,slinhi,slinlogt22,slinhigt22,slinlogt5d22,slinhigt5d22, wncent, wscent, $
	file=datapath+'smart_cycle_edge_linfits_20110330.sav'

;---PLOT 1. - plot upper lower and centroid edges of flux butterfly diagram, also plot detection dots in latitude w/ color for flux
;CROP THIS PLOT:
tcrop=anytim('01-jan-2009')
wcrop=(where(abs(xarr-tcrop) eq min(abs(xarr-tcrop))))[0]
setplotenv,xs=13,ys=13,/ps,file=plotpath+'plot_1_analyze_fluxbflydiag.eps'
setcolors,/sys,/silent,/quiet

;Plot centroid lines
utplot,xarr-mintimseries,nfluxcent,mintimseries,yran=[-60,60],/xsty,xran=[0,tcrop-mintimseries],pos=[.12,.525,.95,.95],xtickname=strarr(10)+' ',xtit='',/nodata,ytit='Latitude',thick=5
oplot,xarr-mintimseries,nfluxcent,color=!gray,ps=10 ;,thick=3
oplot,xarr-mintimseries,sfluxcent,color=!gray,ps=10 ;,thick=3

;Over plot north edges
oplot,xarr[0:nlim]-mintimseries,nfluxmm[0,0:nlim],ps=1,color=!gray
oplot,xarr-mintimseries,nfluxmm[1,*],ps=1,color=!gray
oplot,xarr[0:nlim]-mintimseries,nfluxmmgt22[0,0:nlim],color=!gray,ps=4
oplot,xarr-mintimseries,nfluxmmgt22[1,*],color=!gray,ps=4
;oplot,xarr[0:nlim]-mintimseries,nfluxmmgt5d22[0,0:nlim],color=!gray,ps=2
;oplot,xarr-mintimseries,nfluxmmgt5d22[1,*],color=!gray,ps=2

;Over plot north edge line fits
oplot,xarr[nmin:ncentlim]-mintimseries,nlincent[1]*(xarr[nmin:ncentlim])/secperyr+nlincent[0]
oplot,xarr[nmin:nlim]-mintimseries,nlinlo[1]*(xarr[nmin:nlim])/secperyr+nlinlo[0]
;oplot,xarr[nmin:*]-mintimseries,nlinhi[1]*(xarr[nmin:*])/secperyr+nlinhi[0]
oplot,xarr[nmin:nlim]-mintimseries,nlinlogt22[1]*(xarr[nmin:nlim])/secperyr+nlinlogt22[0]
;oplot,xarr[nmin:*]-mintimseries,nlinhigt22[1]*(xarr[nmin:*])/secperyr+nlinhigt22[0]
;oplot,xarr[nmin:nlim]-mintimseries,nlinlogt5d22[1]*(xarr[nmin:nlim])/secperyr+nlinlogt5d22[0]
;oplot,xarr[nmin:*]-mintimseries,nlinhigt5d22[1]*(xarr[nmin:*])/secperyr+nlinhigt5d22[0]

;Over plot South edges
oplot,xarr[0:slim]-mintimseries,sfluxmm[1,0:slim]-60.,ps=1,color=!gray
oplot,xarr-mintimseries,sfluxmm[0,*]-60.,ps=1,color=!gray
oplot,xarr[0:slim]-mintimseries,sfluxmmgt22[1,0:slim]-60.,color=!gray,ps=4
oplot,xarr-mintimseries,sfluxmmgt22[0,*]-60.,color=!gray,ps=4
;oplot,xarr[0:slim]-mintimseries,sfluxmmgt5d22[1,0:slim]-60.,color=!gray,ps=10
;oplot,xarr-mintimseries,sfluxmmgt5d22[0,*]-60.,color=!gray,ps=10

;Over plot South edge line fits
oplot,xarr[smin:scentlim]-mintimseries,slincent[1]*(xarr[smin:scentlim])/secperyr+slincent[0]
;oplot,xarr[smin:slim]-mintimseries,(slinlo[1]*(xarr[smin:slim])/secperyr+slinlo[0])
oplot,xarr[smin:slim]-mintimseries,(slinhi[1]*(xarr[smin:slim])/secperyr+slinhi[0])
;oplot,xarr[smin:slim]-mintimseries,(slinlogt22[1]*(xarr[smin:slim])/secperyr+slinlogt22[0])
oplot,xarr[smin:slim]-mintimseries,(slinhigt22[1]*(xarr[smin:slim])/secperyr+slinhigt22[0])
;oplot,xarr[smin:slim]-mintimseries,(slinlogt5d22[1]*(xarr[smin:slim])/secperyr+slinlogt5d22[0])-60.
;oplot,xarr[smin:slim]-mintimseries,(slinhigt5d22[1]*(xarr[smin:slim])/secperyr+slinhigt5d22[0])-60.

vline,xarr[wmiss]-mintimseries,color=!white,thick=16
hline,0,color=!gray,thick=10

legend,[textoidl('\Phi')+' Edge',textoidl('\Phi')+' > 1E22 Mx Edge',textoidl('\Phi')+'Centroid','Linear Fits'],psym=[1,4,0,0],color=[!gray,!gray,!gray,0],/top,/right,chars=1.4;,number=[1,1,.5,.5]

fluxstack=[[-sFLUXcont-sFLUXgt22cont-sFLUXgt5d22cont-sFLUXgt23cont],[-nFLUXcont-nFLUXgt22cont-nFLUXgt5d22cont-nFLUXgt23cont]]
;fluxstack=[[sFLUXcont+sFLUXgt22cont+sFLUXgt5d22cont+sFLUXgt23cont],[nFLUXcont+nFLUXgt22cont+nFLUXgt5d22cont+nFLUXgt23cont]]
fluxstack[where(fluxstack eq fluxstack[0,0])]=max(fluxstack)
plot_image,fluxstack[0:wcrop,*],/noerase,pos=[.12,.1,.95,.525],color=255,/xsty,/nosq
contour,([[sFLUXcontsep],[nFLUXcontsep]])[0:wcrop,*],/over,color=0 ;outline "this" solar cycle.
utplot,xarr-mintimseries,nfluxcent,mintimseries,yran=[-60,60],/xsty,xran=[0,tcrop-mintimseries],thick=3,pos=[.12,.1,.95,.525],/nodata,/noerase,ytit='Latitude',xtit='Year (Begins 1 June 1996)'

;shade out prev SC stuff
;polyfill,[-nlinlo[0]/nlinlo[1]*secperyr,-slinhi[0]/slinhi[1]*secperyr,anytim(mintimseries),anytim(mintimseries)]-mintimseries,[0,0,slinhi[1]*(anytim(mintimseries)/secperyr)+slinhi[0],nlinlo[1]*(anytim(mintimseries)/secperyr)+nlinlo[0]], $
;	/data,/LINE_FILL,color=100,spacing=.2

vline,xarr[wmiss]-mintimseries,color=!white,thick=16
hline,0,color=!gray,thick=10
plot_image,rebin([[3,2,1,0],[3,2,1,0]],200,100,/samp),position=[.6,.485,.85,.505],/noerase,yticklen=.0001,xticklen=1,xtickname=strarr(10)+' ',ytickname=strarr(10)+' ',color=0;,/nosq
;plot_image,rebin([[0,1,2,3],[0,1,2,3]],200,100,/samp),position=[.6,.485,.85,.505],/noerase,yticklen=.0001,xticklen=1,xtickname=strarr(10)+' ',ytickname=strarr(10)+' ',color=0;,/nosq
!p.charsize=1.4
xyouts,.6,.458,'   >0         >1E22   >5E22   >1E23 [Mx]',/norm,charsize=1.4

closeplotenv

stop

;Create LATEX table for line fits
nconstrain=2. ;number of fitting contraints [b,m] (?)
slopes=strtrim(string([nlincent[1], nlinlo[1], nlinlogt22[1], slincent[1], slinhi[1], slinhigt22[1]],form='(F7.2)'),2)
yints=strtrim(string([nlincent[0], nlinlo[0], nlinlogt22[0], slincent[0], slinhi[0], slinhigt22[0]],form='(F7.2)'),2)
;include uncertainties
slopes=slopes+' +-'+strtrim(string([nlincenterr[0],nlinloerr[0],0,slincenterr[0],slinhierr[0],0],form='(F7.2)'),2)
yints=yints+' +-'+strtrim(string([nlincenterr[1],nlinloerr[1],0,slincenterr[1],slinhierr[1],0],form='(F7.2)'),2)
chisqred=strtrim(string([ncentchi/(n_elements(wncent)-nconstrain), nlochi/(n_elements(wnhi)-nconstrain), scentchi/(n_elements(wscent)-nconstrain), shichi/(n_elements(wshi)-nconstrain)],form='(F7.2)'),2)
;rowname=['N. Centroid','N. Equatorward','N. Eq. > 1E22 Mx','S. Centroid','S. Equatorward','S. Eq. > 1E22 Mx']
rowname=['','Fit',' $ mX+b $ ','reduced $ \Chi^2 $ ']
;table1data=transpose([[rowname],[' $ '+slopes+'X'+'+'+yints+' $ '],[chisqred]])
table1data=[['','N','S','N','S'],[rowname[1],' $ '+slopes[[0,3,1,4]]+'X'+'+'+yints[[0,3,1,4]]+' $ '],[rowname[2],chisqred]]
;table1cols=['Fit',' $ mX+b $ ','reduced $ \Chi^2 $ ']
table1cols=['','\multicolumn{2}{l|}{Equatorward}','\multicolumn{2}{l|}{Centroid}']
array2textable, table1data, table1cols, 'Linear fits to poleward edges and latitude centroids of hemispheric flux in MF detections.', label='table_1_line_fits', outfile=datapath+'figure_1_linefit_table.txt',/quiet
;---END PLOT 1.--------------------------------------------------------------->


;---PLOT 1A. - Calculate the time of AR emergence for each latitude bin.
xbin=(xarr[1]-xarr[0])/3600./24.
print,anytim(xarr[0],/vms)
x0=xarr[0]
tmin=1979.

;get rid of points wih no edge detected
nsmin=[40,60] ;north then south for minimum index to include in the fit
;if (where(nfluxmm eq missing))[0] ne -1 then nfluxmm[where(nfluxmm eq missing)]=0./0.
;if (where(sfluxmm eq missing))[0] ne -1 then sfluxmm[where(sfluxmm eq missing)]=0./0.
wslo=where(finite(sfluxmm[0,nsmin[1]:*]) eq 1)
wnhi=where(finite(nfluxmm[1,nsmin[0]:*]) eq 1)

;re-do poleward line fits
slinlo=linfit((xarr[nsmin[1]:*])[wslo]/secperyr,(sfluxmm[0,nsmin[1]:*])[wslo]-60.);, MEASURE=fltarr(n_elements(wslo))+.5, chisq=slochi)
nlinhi=linfit((xarr[nsmin[0]:*])[wnhi]/secperyr,(nfluxmm[1,nsmin[0]:*])[wnhi]);, MEASURE=fltarr(n_elements(wnhi))+.5);, chisq=nhichi)

;create extrapolated lines to match with solar cycle contour
ynlinlo=(xarr/secperyr)*nlinlo[1]+nlinlo[0]
ynlinhi=(xarr/secperyr)*nlinhi[1]+nlinhi[0]
yslinlo=(xarr/secperyr)*slinlo[1]+slinlo[0]
yslinhi=(xarr/secperyr)*slinhi[1]+slinhi[0]

nx=n_elements(xarr)
ny=n_elements(yarr)/2.

xind=findgen(nx)

;make test plot
plot_image,[[sfluxcont+sfluxcontsep],[nfluxcont+nfluxcontsep]]
oplot,xind,ynlinlo+60.,color=!red,ps=-4
oplot,xind,ynlinhi+60.,color=!red,ps=-4
oplot,xind,yslinlo+60.,color=!red,ps=-4
oplot,xind,yslinhi+60.,color=!red,ps=-4

ndetmask=fltarr(nx,ny)+1.
sdetmask=fltarr(nx,ny)+1.

spawn,'echo '''' >> '+datapath+'test_mask_detections_range.txt'
for i=0,nx-1 do begin
	if ynlinlo[i] ge 0 then ndetmask[i,0:ynlinlo[i]-1]=0
	if yslinlo[i] le 0 then sdetmask[i,0:yslinlo[i]-1+60.]=0
	if ynlinhi[i] ge 0 then ndetmask[i,ynlinhi[i]:*]=0
	if yslinhi[i] le 0 then sdetmask[i,yslinhi[i]+60.:*]=0
plot_image,[[sdetmask],[ndetmask]]
spawn,'echo ''0:'+strtrim(ynlinlo[i]-1,2)+' , 0:'+strtrim(yslinlo[i]-1+60.,2)+' , '+strtrim(ynlinhi[i],2)+':* , '+strtrim(yslinhi[i]+60.,2)+':*'' >> '+datapath+'test_mask_detections_range.txt'
endfor
;for i=0,nx-1 do begin & if ynlinlo[i] ge 0 then ndetmask[i,0:ynlinlo[i]-1]=0 & if yslinlo[i] le 0 then sdetmask[i,0:yslinlo[i]-1+60.]=0 & if ynlinhi[i] ge 0 then ndetmask[i,ynlinhi[i]:*]=0 & if yslinhi[i] le 0 then sdetmask[i,yslinhi[i]+60.:*]=0 & plot_image,[[sdetmask],[ndetmask]] & spawn,'echo ''0:'+strtrim(ynlinlo[i]-1,2)+' , 0:'+strtrim(yslinlo[i]-1+60.,2)+' , '+strtrim(ynlinhi[i],2)+':* , '+strtrim(yslinhi[i]+60.,2)+':*'' >> '+datapath+'test_mask_detections_range.txt' & endfor

maskcontigblob=[[sfluxcont],[nfluxcont]]
maskcontigblob[where([[sdetmask],[ndetmask]] eq 1)]=1
;masksepblob = LABEL_REGION(maskcontigblob) ;
masksepblob = LABEL_REGION([[fltarr(nx+2)],[([transpose(fltarr(ny*2.)),maskcontigblob,transpose(fltarr(ny*2.))])],[fltarr(nx+2)]])
masksepblob=masksepblob[1:nx,1:ny*2.]

;Make time of emergence plot
timearcycle=[[sfluxcont],[nfluxcont]]
timearcycle[where(masksepblob gt 1)]=0.

;find time difference between minimum and maximum emergence
dtemergelat=fltarr(ny*2)
for i=0,ny*2-1 do begin
	thissten=minmax(where(timearcycle[*,i] eq 1))
	if total(thissten) eq 0 then dtemergelat[i]=0 $
		else dtemergelat[i]=(thissten[1]-thissten[0]+1.)*(xarr[1]-xarr[0])/3600./24./365.
endfor
;for i=0,ny*2-1 do begin & thissten=minmax(where(timearcycle[*,i] eq 1)) & if total(thissten) eq 0 then dtemergelat[i]=0 else dtemergelat[i]=(thissten[1]-thissten[0]+1.)*(xarr[1]-xarr[0])/3600./24./365. & endfor

setplotenv,/ps,file=plotpath+'plot_1a_analyze_timearemerge.eps'
setcolors,/sile,/sys
;estimate flux at each latitude
plot,yarr,total([[sFLUXIMAGE],[nFLUXIMAGE]]*1d16*timearcycle,1), /xsty, xran=[-60,60], $
	xticklen=.0001, yticklen=.0001, xtickname=strarr(10)+' ', ytickname=strarr(10)+' ', /ysty, yran=[0,15d23], color=!gray, thick=5, ps=10, $
	pos=[.1,.15,.85,.95]
plot,yarr,total(timearcycle,1)*(xarr[1]-xarr[0])/3600./24./365.,ps=10,yran=[0,15], $
	ytitle='MF Emergence [Years]', xtitle='Latitude',thick=10, /xsty, xran=[-60,60],/noerase, $
	pos=[.1,.15,.85,.95]
oplot,yarr,dtemergelat,ps=10
vline,0,lines=2,yran=[0,11.5]
;oplot,yarr,total([[sdetmask],[ndetmask]],1)*(xarr[1]-xarr[0])/3600./24./365.,ps=10,color=150
axis,yaxis=1, /ysty, yran=[0,15d23],color=!gray,ytit=textoidl('\Phi_{TOT}')+' [Mx]'
legend, /left,/top,[textoidl('\Sigma t')+' Emergence',textoidl('\Delta t')+' Start-End',textoidl('\Sigma\Phi_{TOT}')],lines=[0,0,0],color=[0,0,!gray],thick=[5,10,5]
closeplotenv


stop

;---END PLOT 1A.-------------------------------------------------------------->


;PLOT 2 JUMP POINT
gotoplot2:


;CREATE PLOT OF TIME EACH LAT PRODUCES EMERGENT FLUX-------------------------->

;DO THIS!!!
;nFLUXcontsep=smart_cont_sep(nFLUXcont, contlevel=.5, vthresh=1, areathresh=1)
;sFLUXcontsep=smart_cont_sep(sFLUXcont, contlevel=.5, vthresh=1, areathresh=1)

;plot_image,nFLUXcontsep







stop

;CREATE YEAR BIN TSERIES------------------------------------------------------>

;Define epochs of cycle	
trise=anytim('1-jan-1997')
tmax=anytim('1-jan-2000')
tdecay=anytim('1-jan-2003')
tend=anytim('1-jan-2009')

;Create year binned time series
if not keyword_set(res4store) then begin
	tstart=anytim('1-jan-1996')
	smart_bin_time, arstr_arr.area, ARSTR_ARR_tim, areayrbin, tareayrbin,/year,/total, tstart=tstart
	smart_bin_time, arstr_arr.bflux, ARSTR_ARR_tim, fluxyrbin, tfluxyrbin,/year,/total, tstart=tstart
	smart_bin_time, arstr_arr.area, ARSTR_ARR_tim, areameanyrbin, tareayrbin,/year,/mean, tstart=tstart
	smart_bin_time, arstr_arr.bflux, ARSTR_ARR_tim, fluxmeanyrbin, tfluxyrbin,/year,/mean, tstart=tstart
	
	save,areayrbin,tareayrbin,fluxyrbin,tfluxyrbin,areameanyrbin,fluxmeanyrbin,file=datapath+'smart_yearbins_20110330.sav'
endif else restore,datapath+'smart_yearbins_20110330.sav',/ver


;Create AR distributions for each epoch
wrise=where(ARSTR_ARR_tim ge trise and ARSTR_ARR_tim lt tmax)

;Filter out old cycle regions
restore,datapath+'smart_cycle_edge_linfits_20110330.sav',/ver
wnewcyc=[where(arstr_arr[wrise].hglat ge 0 and arstr_arr[wrise].hglat gt (nlinlo[1]*ARSTR_ARR_tim[wrise]/secperyr+nlinlo[0])), $
		where(arstr_arr[wrise].hglat lt 0 and arstr_arr[wrise].hglat lt (slinhi[1]*ARSTR_ARR_tim[wrise]/secperyr+slinhi[0]))]
wrise=wrise[wnewcyc]

wmax=where(ARSTR_ARR_tim ge tmax and ARSTR_ARR_tim lt tdecay)
wdecay=where(ARSTR_ARR_tim ge tdecay and ARSTR_ARR_tim lt tend)

ardistrise=arstr_arr[wrise]
ardistmax=arstr_arr[wmax]
ardistdecay=arstr_arr[wdecay]







;---PLOT 2. - separate distributions of flux and area for 3 phases of solar cycle

;Shift plots by half a bin so the bars cover the full bin instead of being centered on them
halfbin=.5*365.*3600.*24. ;seconds per half year
areameanyrbin=[areameanyrbin[0],areameanyrbin]
tareayrbin=[(tareayrbin)[0],tareayrbin+halfbin]
fluxmeanyrbin=[fluxmeanyrbin[0],fluxmeanyrbin]
tfluxyrbin=[(tfluxyrbin)[0],tfluxyrbin+halfbin]

reftim=min(tareayrbin)

;USE HISTOGRAM FOR DISTRIBUTION PLOTS
;AREA:
ardistareamax = HISTOGRAM( ardistmax.area, BINSIZE=2500., locations=xareamax)
ardistarearise = HISTOGRAM( ardistrise.area, BINSIZE=2500.,locations=xarearise)
ardistareadecay = HISTOGRAM( ardistdecay.area, BINSIZE=2500.,locations=xareadecay)

;FLUX:
ardistfluxmax = HISTOGRAM( ardistmax.bflux*1d16, BINSIZE=3d21, locations=xfluxmax)
ardistfluxrise = HISTOGRAM( ardistrise.bflux*1d16, BINSIZE=3d21,locations=xfluxrise)
ardistfluxdecay = HISTOGRAM( ardistdecay.bflux*1d16, BINSIZE=3d21,locations=xfluxdecay)

;Do power-law fitting
;AREA:
areafitlim=[1d3,3d4]
plawareamax=fit_plaw(xareamax, ardistareamax, yfit=fitareamax, xfit=xfitareamax, perror=fitareamaxerr, limxfit=areafitlim, chisqred=chiareamax)
plawarearise=fit_plaw(xarearise, ardistarearise, yfit=fitarearise, xfit=xfitarearise, perror=fitareariseerr, limxfit=areafitlim, chisqred=chiarearise)
plawareadecay=fit_plaw(xareadecay, ardistareadecay, yfit=fitareadecay, xfit=xfitareadecay, perror=fitareadecayerr, limxfit=areafitlim, chisqred=chiareadecay)

;FLUX:
fluxfitlim=[1d21,3.5d22]
plawfluxmax=fit_plaw(xfluxmax, ardistfluxmax, yfit=fitfluxmax, xfit=xfitfluxmax, perror=fitfluxmaxerr, limxfit=fluxfitlim, chisqred=chifluxmax)
plawfluxrise=fit_plaw(xfluxrise, ardistfluxrise, yfit=fitfluxrise, xfit=xfitfluxrise, perror=fitfluxriseerr, limxfit=fluxfitlim, chisqred=chifluxrise)
plawfluxdecay=fit_plaw(xfluxdecay, ardistfluxdecay, yfit=fitfluxdecay, xfit=xfitfluxdecay, perror=fitfluxdecayerr, limxfit=fluxfitlim, chisqred=chifluxdecay)

;START PLOTTING--------------->
savpos=!p.position
setplotenv,xs=15,ys=18,/ps,file=plotpath+'plot_2_phase_cycle_ar_prop.eps'
setcolors,/sys,/silent,/quiet

;Yearly averaged area and flux with vertical lines showing breaks between phases in cycle
utplot,tareayrbin-reftim,areameanyrbin,reftim,ps=10,position=[.12,.83,.97,.99],xtit='',xtickname=strarr(10)+' ',ytitle='<Area> [Mm'+textoidl('^2')+']',/xsty;ytickformat='exponent',ytickformat='(E9.1)'
vline,[trise,tmax,tdecay,tend]-reftim,thick=3,color=!gray

;indicate peaks on year bin plots
plotsym,2,5,thick=5
plots,[anytim('1-jun-2000'),anytim('1-jun-2002'),anytim('1-jun-2004')]-reftim,[6000,6000,6000],ps=8

xyouts,anytim('1-jun-1997')-reftim,1000,'Rise',/data
xyouts,anytim('1-jun-2000')-reftim,1000,'Plateau',/data
xyouts,anytim('1-jun-2004')-reftim,1000,'Decay',/data

utplot,tfluxyrbin-reftim,fluxmeanyrbin*1d16,reftim,ps=10,position=[.12,.67,.97,.83],/noerase,ytitle='<'+textoidl('\Phi_{TOT}')+'> [Mx]',/xsty,xtit='Year (Begins 1 January 1996)'
vline,[trise,tmax,tdecay,tend]-reftim,thick=3,color=!gray

;indicate peaks on year bin plots
plots,[anytim('1-jun-1999'),anytim('1-jun-2002'),anytim('1-jun-2004')]-reftim,[5d21,5d21,5d21],ps=8,thick=2

;plot area distributions here
!p.position=[.12,.37,.97,.62]
plot,xareamax,ardistareamax,ps=10,/xlog,/ylog,/xsty,/ysty,xran=[5d2,3d5],yran=[1,1d6],ytitle='# Detections',chars=!p.charsize,xtit='Area [Mm'+textoidl('^2')+']',/noerase
oplot,xareamax,ardistareamax,ps=4
oplot,xarearise,ardistarearise,ps=10,lines=2
oplot,xareadecay,ardistareadecay,ps=10,color=!gray

;plot area P-Law fits
oplot,xfitareamax, fitareamax,color=!red,ps=-1,thick=5
;oplot,xfitarearise, fitarearise,color=!red,ps=-1
;oplot,xfitareadecay, fitareadecay,color=!red,ps=-1

;bin=.02
;plot_hist,alog10(ardistmax.area),ytitle='# Detections',/noerase,bin=bin,/log,xran=[3.5,5.5],/xsty,chars=!p.charsize,xtit='Log'+textoidl('_{10}(Area)')+' [Mm'+textoidl('^2')+']'
;plot_hist,alog10(ardistrise.area),/oplot,bin=bin,/log,lines=2
;plot_hist,alog10(ardistdecay.area),/oplot,bin=bin,/log,color=!gray

legend,['Rise','Plateau','Decay'],lines=[2,0,0],color=[!black,!black,!gray],/right,/top

;plot flux distributions here
!p.position=[.12,.06,.97,.31]
plot,xfluxmax,ardistfluxmax,ps=10,/xlog,/ylog,/xsty,/ysty,xran=[1d20,5d23],yran=[1,1d6],ytitle='# Detections',chars=!p.charsize,xtit=textoidl('\Phi_{TOT}')+' [Mx]',/noerase
oplot,xfluxmax,ardistfluxmax,ps=4
oplot,xfluxrise,ardistfluxrise,ps=10,lines=2
oplot,xfluxdecay,ardistfluxdecay,ps=10,color=!gray

;plot flux P-Law fits
oplot,xfitfluxmax, fitfluxmax,color=!red,ps=-1,thick=5
;oplot,xfitfluxrise, fitfluxrise,color=!red,ps=-1
;oplot,xfitfluxdecay, fitfluxdecay,color=!red,ps=-1

;plot_hist,alog10(ardistmax.bflux*1d16),ytitle='# Detections',/noerase,bin=.01,/log,xran=[21,23.5],/xsty,chars=!p.charsize,xtit='Log'+textoidl('_{10}(\Phi_{TOT})')+' [Mx]'
;plot_hist,alog10(ardistrise.bflux*1d16),/oplot,bin=bin,lines=2
;plot_hist,alog10(ardistdecay.bflux*1d16),/oplot,bin=bin,color=!gray

;Create LATEX table for power-law fits
nconstrain=2. ;number of fitting contraints [b,m] (?)
rowname=['Rise','Plateau','Decay']
areaplawcol='$ '+strtrim(string([plawarearise[1],plawareamax[1],plawareadecay[1]],form='(F7.2)'),2)+' $'
areaplawcolerr='$ '+strtrim(string([fitareariseerr[1],fitareamaxerr[1],fitareadecayerr[1]],form='(F7.2)'),2)+' $'
fluxplawcol='$ '+strtrim(string([plawfluxrise[1],plawfluxmax[1],plawfluxdecay[1]],form='(F7.2)'),2)+' $'
fluxplawcolerr='$ '+strtrim(string([fitfluxriseerr[1],fitfluxmaxerr[1],fitfluxdecayerr[1]],form='(F7.2)'),2)+' $'
areachisqcol='$ '+strtrim(string([chiarearise,chiareamax,chiareadecay],form='(F7.2)'),2)+' $'
fluxchisqcol='$ '+strtrim(string([chifluxrise,chifluxmax,chifluxdecay],form='(F7.2)'),2)+' $'
table1data=transpose([[rowname],[areaplawcol+' $ \pm $ '+areaplawcolerr],[areachisqcol],[fluxplawcol+' $ \pm $ '+fluxplawcolerr],[fluxchisqcol]])
table1cols=['Cycle Phase','Area',' $ \Chi^2 $ ','Flux',' $ \Chi^2 $ ']
array2textable, table1data, table1cols, 'Power-law fits to distributions of MF area and flux.', label='table_2_line_fits', outfile=datapath+'figure_3_plawfit_table.txt',/quiet


!p.position=savpos

closeplotenv

stop

;---END PLOT 2.--------------------------------------------------------------->


;PLOT 3. - Combine PFSS/LASCO comparisons into a mosaic

;See ~/science/procedures/pfss/compare_pfss_lasco.pro

;---END PLOT 3.--------------------------------------------------------------->


;ANALYZE MAGNETIC BUTTERFLY DIABGRAM------------------------------------------>
restore,datapath+'smart_cycle_edge_linfits_20110330.sav' ;to get WNCENT and WSCENT
restore,'~/science/data/cycle23_sav/carrington_mdi/cycle23_butterfly_flux_maps_final.sav',/ver
restore,'~/science/papers/active_regions_2_cycle/data/butterfly_diags_lat_27day_bin.sav',/ver
restore,datapath+'pfss_monthly_phiab_plot.sav',/ver
;rebin butterfly diagram to 27 day bins by 1 degree
bin=27. ;in days
magrebin=rebin1d(MEANNETB60,bin,direct=1)
wbad=where(total(magrebin[*,24:154],2) gt 200 or total(magrebin[*,24:154],2) lt (-200))
magrebin[wbad,*]=0./0.
timrebin=findgen(floor(n_elements(TIMARR)/bin))*bin*3600.*24.+min(TIMARR) 

;north pole average field
npole6065=average(magrebin[*,150:154],2)
npole6570=average(magrebin[*,155:159],2)
npole7075=average(magrebin[*,160:164],2)

;south pole average field
spole6065=average(magrebin[*,25:29],2)
spole6570=average(magrebin[*,20:24],2)
spole7075=average(magrebin[*,15:19],2)

;mask flux signed butterfly diagram to pull out N and S high and low-latitude flux imbalances 
nbins=n_elements(xarr)
nfluxmask=fltarr(nbins)
nfluxmask[wncent]=1
sfluxmask=fltarr(nbins)
sfluxmask[wscent]=1
wncentmask=where(nfluxmask eq 0)
wscentmask=where(sfluxmask eq 0)
FLUXSIGNEDhigh=FLUXSIGNEDIMAGE
FLUXSIGNEDhigh[wncentmask,60:*]=0./0.
FLUXSIGNEDhigh[wscentmask,0:59]=0./0.
FLUXSIGNEDlow=FLUXSIGNEDhigh

for i=0,nbins-1 do begin
	if nfluxmask[i] then begin
		FLUXSIGNEDhigh[i,60:60.+nfluxcent[i]-1.]=0.
		FLUXSIGNEDlow[i,60.+nfluxcent[i]:*]=0.
	endif
	if sfluxmask[i] then begin
		FLUXSIGNEDhigh[i,60.+sfluxcent[i]:59]=0.
		FLUXSIGNEDlow[i,0:60.+sfluxcent[i]-1.]=0.
	endif
endfor

;High-latitude AR net flux
nfluxnethigh=total(FLUXSIGNEDhigh[*,60:*],2)
sfluxnethigh=total(FLUXSIGNEDhigh[*,0:59],2)

;Low-latitude AR net flux
nfluxnetlow=total(FLUXSIGNEDlow[*,60:*],2)
sfluxnetlow=total(FLUXSIGNEDlow[*,0:59],2)

;Find where times best match for two images (one starts earlier, one ends later)
wcliptim0=(where(abs(timrebin-min(xarr)) eq min(abs(timrebin-min(xarr)))))[0]
wcliptim1=(where(abs(xarr-max(timrebin)) eq min(abs(xarr-max(timrebin)))))[0]
magrebincrop=magrebin[wcliptim0:*,*]
timrebincrop=timrebin[wcliptim0:*]
FLUXSIGNEDcrop=FLUXSIGNEDIMAGE[0:wcliptim1,*]
xarrcrop=xarr[0:wcliptim1]

nbinscrop=n_elements(xarrcrop)

;Pull out mean B field north and south in unipolar flows
nbsignedunipole=fltarr(nbinscrop)+0./0.
sbsignedunipole=fltarr(nbinscrop)+0./0.
nbmask=fltarr(nbinscrop)
nbmask[wnhi]=1
sbmask=fltarr(nbinscrop)
sbmask[wslo]=1

;Crop to match size of B signed array
nbmask=nbmask[0:wcliptim1] & sbmask=sbmask[0:wcliptim1]
nfluxmmcrop=nfluxmm[*,0:wcliptim1] & sfluxmmcrop=sfluxmm[*,0:wcliptim1]

wnhimask=where(nbmask eq 0)
wslomask=where(sbmask eq 0)
magrebinhigh=magrebincrop
magrebinhigh[wnhimask,90:*]=0./0.
magrebinhigh[wslomask,0:89]=0./0.
magrebinhighmask=fltarr((size(magrebinhigh))[1],(size(magrebinhigh))[2])

for i=0,nbinscrop-1 do begin
	if nbmask[i] then begin
		nbsignedunipole[i]=mean(magrebinhigh[i,(nfluxmmcrop[1,i]+90.)>0:(nfluxmmcrop[1,i]+90.+10.)>0]) ;added ">0" because nfluxmmcrop has -NaN values! Need to check this...
		magrebinhighmask[i,(nfluxmmcrop[1,i]+90.)>0:(nfluxmmcrop[1,i]+90.+10.)>0]=1
	endif
	if sbmask[i] then begin
		sbsignedunipole[i]=mean(magrebinhigh[i,(sfluxmmcrop[0,i]+30.-10.)>0:(sfluxmmcrop[0,i]+30.)>0])
		magrebinhighmask[i,(sfluxmmcrop[0,i]+30.-10.)>0:(sfluxmmcrop[0,i]+30.)>0]=1
	endif
endfor

;plot_image,(magrebinhigh*magrebinhighmask)<10>(-10),/nosq
;stop

;---PLOT 4. - Butterfly avg signed B with flux imbalance butterfly------------>

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_4_magbutt_arsignflux.eps'

loadct,0
mintim=min(timrebin)

;plot b field butterfly diagram
!p.position=[.12,.525,.95,.95]
;plot_image,magrebin[*,24:153]<3>(-3),/nosq,color=255,xticklen=.0001,yticklen=.0001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' '
;utplot,timrebin-mintim,findgen(n_elements(timrebin)),mintim,/nodat,/noerase,/xsty,yran=[-90+24,89-24.],/ysty,ytit='Latitude'
plot_image,magrebincrop[*,24:153]<3>(-3),/nosq,color=255,xticklen=.0001,yticklen=.0001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' '
utplot,timrebincrop-mintim,findgen(n_elements(timrebincrop)),mintim,/nodat,/noerase,/xsty,yran=[-90+24,89-24.],/ysty,ytit='Latitude'
xyouts,.7,.87,'<B'+textoidl('_{SIGNED}')+'> [G]',/norm,color=255
colorbar,range=[-3,3],pos=[.65,.92,.90,.95],color=255,chars=!p.charsize ;,/invert

;plot AR net flux
!p.position=[.12,.1,.95,.525]
plot_image,FLUXSIGNEDcrop<150000>(-150000),/nosq,/noerase,xticklen=.0001,yticklen=.0001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' ',ytit='Latitude'
utplot,XARRcrop-mintim,findgen(n_elements(timrebin)),mintim,/nodat,/noerase,/xsty,yran=minmax(yarr),/ysty
xyouts,.66,.445,'Detection '+textoidl('\Phi_{NET}')+' 10!u20!n [MX]',/norm,color=0
colorbar,range=[-15,15], DIVISIONS=0.5, $
	pos=[.65,.495,.90,.525],color=0,chars=!p.charsize ;,/invert

closeplotenv

;---END PLOT 4.--------------------------------------------------------------->


;PLOT 5. - Zonal polarity imbalances and netflux butterfly-------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_5_fluximbal_hi_lo.eps'
setcolors,/sys,/silent,/quiet
mintim=min(xarr)
xran=[anytim('1-jan-1997'),anytim('1-jan-2006')]-mintim

;Plot with hi and low lat bounds of AR butterfly diagram
utplot,xarr-mintim,nfluxnethigh*1d16,mintim,ps=10,ytit=textoidl('\Sigma\Phi_{NET}')+' [Mx] Mid Lat.',pos=[.12,.666,.95,.95],yran=[-5d22,5d22],/ysty,xtit='',xtickname=strarr(10)+' ',xran=xran,/xsty
oplot,xarr-mintim,nfluxnethigh*1d16,color=!red,ps=10
oplot,xarr-mintim,sfluxnethigh*1d16,color=!blue,ps=10
hline,0,col=!gray
utplot,xarr-mintim,nfluxnetlow*1d16,mintim,ps=10,ytit=textoidl('\Sigma\Phi_{NET}')+' [Mx] Low Lat.',pos=[.12,.383,.95,.666],/noerase,yran=[-5d22,5d22],/ysty,xtit='',xtickname=strarr(10)+' ',xran=xran,/xsty
oplot,xarr-mintim,nfluxnetlow*1d16,color=!red,ps=10
oplot,xarr-mintim,sfluxnetlow*1d16,color=!blue,ps=10
hline,0,col=!gray

;Average B signed in some latitude range in flux butterfly diagram
utplot,timrebincrop-mintim,nbsignedunipole,mintim,ps=10,pos=[.12,.1,.95,.383],ytit='<B'+textoidl('_{SIGNED}')+'> High Lat.',/noerase,yran=[-10,10],/ysty,xran=xran,/xsty
oplot,timrebincrop-mintim,nbsignedunipole,color=!red,ps=10
oplot,timrebincrop-mintim,sbsignedunipole,color=!blue,ps=10
hline,0,col=!gray

closeplotenv

stop


;---END PLOT 5.--------------------------------------------------------------->


;PLOT 6. - plot polar fields and dipole sphere harmonics---------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_6_polarb_sphereharm.eps'
setcolors,/sys,/silent,/quiet
mintim=min(tflux27dseries)

linecol=0

;polar plot
setcolors,/sys
utplot,timrebin-mintim,npole6065,mintim,yran=[-10,10], pos=[.12,.525,.95,.95],ps=10,xran=minmax(timrebin)-mintim,/xsty,ytit='Polar Fields [G]',xtit='',xtickname=strarr(10)+' '
oplot,timrebin-mintim,npole6065,color=!red,ps=10
oplot,timrebin-mintim,npole6570,color=!red,ps=10,lines=2
;oplot,timrebin-mintim,npole7075,color=!red,ps=10
oplot,timrebin-mintim,spole6065,color=!blue,ps=10
oplot,timrebin-mintim,spole6570,color=!blue,ps=10,lines=2
;oplot,timrebin-mintim,spole7075,color=!blue,ps=10
hline,0,color=!gray

legend,['Lat. 60-65','Lat. 65-70'],/top,/right,color=[linecol,linecol],lines=[0,2]


;spherical harmonic plot
;pure dipole
utplot,ABTIMARR-mintim,reform(COEFF01),mintim,pos=[.12,.1,.95,.525],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-2,2],ytit='Harmonic Coefficients'
;two bands in each hemisphere
oplot,ABTIMARR-mintim,reform(COEFF03),ps=10,color=!gray
;three bands in each hemisphere (AR, decay, pole)
oplot,ABTIMARR-mintim,reform(COEFF05),ps=10,lines=2

;oplot,ABTIMARR-mintim,reform(COEFF02),ps=10,color=!red ;non hemisherically symmetric
;oplot,ABTIMARR-mintim,reform(COEFF04),ps=10,color=!green ;non hemisherically symmetric
hline,0,/gray

legend,[textoidl('C_{1,0}'),textoidl('C_{3,0}'),textoidl('C_{5,0}')],/bottom,/right,color=[linecol,!gray,linecol],lines=[0,0,2]

closeplotenv

;See sphere_harminic.pro for harmonic plots

;See carrington_mdi_analyze.pro for polar plots

stop

;---END PLOT 6.--------------------------------------------------------------->


;PLOT 2.5 - AR complexes

plot_hist,alog10(arstr_arr[wcomplex].bflux*1d16),bin=.1
wcomplex=where(arstr_arr.extstr.HGLONWIDTH gt 50)      
plot_hist,alog10(arstr_arr[wcomplex].bflux*1d16),bin=.1,/oplot,lines=2
wcomplex=where(arstr_arr.extstr.HGLONWIDTH gt 60)                     
plot_hist,alog10(arstr_arr[wcomplex].bflux*1d16),bin=.1,/oplot,lines=3
wcomplex=where(arstr_arr.extstr.HGLONWIDTH gt 30)                     
plot_hist,alog10(arstr_arr[wcomplex].bflux*1d16),bin=.1,/oplot,lines=4
wcomplex=where(arstr_arr.extstr.HGLONWIDTH gt 20)                     
plot_hist,alog10(arstr_arr[wcomplex].bflux*1d16),bin=.1,/oplot,lines=2

;wfluxgt=where(alog10(arstr_arr.bflux*1d16) gt 23.5)
;print,arstr_arr[wfluxgt].time
; 4-Aug-1996 01:39:04.602 11-Aug-1996 04:51:04.479 29-Aug-1996 03:15:04.248
;31-Aug-1996 04:51:04.239  1-Sep-1996 04:51:04.236 22-Sep-1996 06:27:04.219
;25-Sep-1996 04:51:04.221 27-Sep-1996 04:51:04.228  5-Feb-1999 17:39:02.128
; 5-Feb-1999 17:39:02.128 24-Feb-1999 22:27:02.148 24-Feb-1999 22:27:02.148
;12-Mar-1999 17:39:02.304 30-Oct-2002 03:15:01.376 30-Dec-2004 06:27:02.109
;30-Dec-2004 06:27:02.109 13-Mar-2005 01:10:03.532 13-Mar-2005 01:10:03.532
;13-Mar-2005 01:10:03.532 13-Mar-2005 01:10:03.532
;print,arstr_arr[wfluxgt].extstr.HGLONWIDTH
;       29.347948       12.386070       28.941887       30.572639
;       33.824963       13.708313       34.178452       38.541359
;       20.401094       19.064800       16.719902       16.009506
;       43.188049       64.828934       61.405422       17.316025
;       19.812109       20.218603       22.092505       18.057243

;plot flux against area
;maximize flux, psl length, minimize area ->most complex, compact region


;PLOT 2.6 - Latitude oscillation

wargt1d21_5=where(arstr_arr.bflux*1d16 gt 10.^(21.5) and arstr_arr.hglat gt 0)
wargt1d22=where(arstr_arr.bflux*1d16 gt 10.^(22) and arstr_arr.hglat gt 0)
wargt1d22_5=where(arstr_arr.bflux*1d16 gt 10.^(22.5) and arstr_arr.hglat gt 0)
wargt1d23=where(arstr_arr.bflux*1d16 gt 10.^(23) and arstr_arr.hglat gt 0)

utplot,xarr-mintim,nfluxcent,mintim,ps=10,xran=[anytim('1-jan-1997'),anytim('1-jan-2006')]-mintim,/xsty
oplot,arstr_arr_tim[wargt1d21_5]-mintim,arstr_arr[wargt1d21_5].hglat,ps=4,color=!blue
oplot,arstr_arr_tim[wargt1d22]-mintim,arstr_arr[wargt1d22].hglat,ps=4,color=!yellow
oplot,arstr_arr_tim[wargt1d22_5]-mintim,arstr_arr[wargt1d22_5].hglat,ps=4,color=!orange
oplot,arstr_arr_tim[wargt1d23]-mintim,arstr_arr[wargt1d23].hglat,ps=4,color=!red
oplot,xarr-mintim,nfluxcent,ps=10,thick=4,color=255
oplot,xarr-mintim,nfluxcent,ps=10,thick=2,color=0

window_capture,file=plotpath+'ar_latitude_oscillation',/png

;PLOT 2.7 - Example plots of various AR fluxes

stop

end