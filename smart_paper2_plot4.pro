;2013-08-07
;Remake plot 4 in the paper

pro smart_paper2_plot4

datapath='./data/'
plotpath='./plots/20110330/'

restore,datapath+'smart_cycle_edge_linfits_20110330.sav',/ver ;to get WNCENT and WSCENT
restore,datapath+'cycle23_butterfly_flux_maps_final.sav',/ver
restore,datapath+'butterfly_diags_lat_27day_bin.sav',/ver
restore,datapath+'pfss_monthly_phiab_plot.sav',/ver
restore,datapath+'wnhi_wslo_variables',/ver

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
;ARs
mfcyclemask=[[sFLUXCONTSEP],[nFLUXCONTSEP]]
mffluxsignimg=mfcyclemask*FLUXSIGNEDIMAGE
mffluxsignn=total(mffluxsignimg[*,60:*],2)*1d16
mffluxsigns=total(mffluxsignimg[*,0:59],2)*1d16

;determine <B> flux at high and low latitudes
;construct mask that fits over signed B mask
bsigncyclemask=([[fltarr(197,30)],[mfcyclemask],[fltarr(197,30)]])[0:188,*]
imgsz=size(bsigncyclemask,/dim)

;get flux centroid pixel positions corresponding to the bsigned map
bsignnfluxcentpx=round(nfluxcent[0:189])+90
bsignsfluxcentpx=round(sfluxcent[0:189])+90

stop

;extend the array in time to start and end
wncycle=where(bsignnfluxcentpx gt 0 and finite(bsignnfluxcentpx) eq 1)
wminmaxncycle=minmax(wncycle)
bsignnfluxcentpx[0:wminmaxncycle[0]]=bsignnfluxcentpx[wminmaxncycle[0]]
bsignnfluxcentpx[wminmaxncycle[1]:*]=bsignnfluxcentpx[wminmaxncycle[1]]

wscycle=where(bsignsfluxcentpx gt 0 and finite(bsignsfluxcentpx) eq 1)
wminmaxscycle=minmax(wscycle)
bsignsfluxcentpx[0:wminmaxscycle[0]]=bsignsfluxcentpx[wminmaxscycle[0]]
bsignsfluxcentpx[wminmaxscycle[1]:*]=bsignsfluxcentpx[wminmaxscycle[1]]

;invert the bsigned masks to get a mask of the high-latitude pixels
bsignhilatmask=fltarr(imgsz[0],imgsz[1])
bsignhilatmask[where(bsigncyclemask eq 0)]=1.

;get rid of stuff below cycle flux centroid position 
for i=0,imgsz[0]-1 do bsignhilatmask[i,90:bsignnfluxcentpx[i]]=0.
for i=0,imgsz[0]-1 do bsignhilatmask[i,bsignsfluxcentpx[i]:89]=0.

;limit to +- 54 degrees
bsignhilatmask[*,90+54:*]=0
bsignhilatmask[*,0:90-54-1]=0

;create a surface area map that will yield flux in the end
latarr=findgen(imgsz[1])-90.
colatarr=reverse(findgen(imgsz[1]))
areaband=2.*!pi*wcs_rsun(unit='cm')^2.*(cos(colatarr*!dtor)-cos((colatarr+1.)*!dtor))
areabandimg=rebin(transpose(areaband),imgsz[0],imgsz[1])
areabandhi=areabandimg*bsignhilatmask
areabandlo=areabandimg*bsigncyclemask

;Determine low lat bsign flux
bsignfluxnlo=total((areabandlo*magrebin)[*,90:*],2)
bsignfluxslo=total((areabandlo*magrebin)[*,0:89],2)

;Determine high lat bsign flux
bsignfluxnhi=total((areabandhi*magrebin)[*,90:*],2)
bsignfluxshi=total((areabandhi*magrebin)[*,0:89],2)

bsignxarr=xarr[0:188]

;replace positions with missing data with NaNs for low lat imbalances
wnan=where(finite(bsignfluxslo) ne 1)
mffluxsignn[wnan]=1/0.
mffluxsigns[wnan]=1/0.

;PLOT 5. - Zonal polarity imbalances and netflux butterfly-------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_5_fluximbal_hi_lo.eps'
setcolors,/sys,/silent,/quiet
mintim=min(xarr)
xran=[anytim('1-jan-1997'),anytim('1-jan-2010')]-mintim

;Plot with hi and low lat bounds of AR butterfly diagram
utplot,bsignxarr-mintim, bsignfluxnhi/1d22,mintim,ps=10,ytit=textoidl('\Phi')+'!d<B>!n ['+textoidl('\times10^{22}')+' Mx] High Lat.',pos=[.12,.525,.95,.95],yran=[-4,4],/ysty,xtit='',xtickname=strarr(10)+' ',xran=xran,/xsty
oplot,bsignxarr-mintim,bsignfluxnhi/1d22,color=!red,ps=10
oplot,bsignxarr-mintim,bsignfluxshi/1d22,color=!blue,ps=10
hline,0,col=!gray
utplot,xarr-mintim,mffluxsignn/1d22,mintim,ps=10,ytit=textoidl('\Phi')+'!dNET,MF!n ['+textoidl('\times10^{22}')+' Mx] Low Lat.',pos=[.12,.1,.95,.525],/noerase,yran=[-4,3.99],/ysty,xran=xran,/xsty
oplot,xarr-mintim,mffluxsignn/1d22,color=!red,thick=5,ps=10
oplot,xarr-mintim,mffluxsigns/1d22,color=!blue,thick=5,ps=10
oplot,bsignxarr-mintim,bsignfluxnlo/1d22,color=!red,thick=2,lines=2,ps=10
oplot,bsignxarr-mintim,bsignfluxslo/1d22,color=!blue,thick=2,lines=2,ps=10
hline,0,col=!gray

closeplotenv


;---END PLOT 5.--------------------------------------------------------------->





;Remake plot with N on top panel and S in bottom panel------------------------>


setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_5_2_fluximbal_hi_lo.eps'
setcolors,/sys,/silent,/quiet
mintim=min(xarr)
xran=[anytim('1-jan-1997'),anytim('1-jan-2011')]-mintim

;Plot with hi and low lat bounds of AR butterfly diagram
utplot,bsignxarr-mintim, bsignfluxnhi/1d22,mintim,ps=10,ytit='North '+textoidl('\Phi_{NET}')+' ['+textoidl('\times10^{22}')+' Mx]',pos=[.12,.525,.95,.95],yran=[-4,4],/ysty,xtit='',xtickname=strarr(10)+' ',xran=xran,/xsty, thick=10
oplot,bsignxarr-mintim,bsignfluxnhi/1d22,color=!red,ps=10, thick=10

oplot,xarr-mintim,mffluxsignn/1d22,color=!red,thick=2,ps=10

oplot,bsignxarr-mintim,bsignfluxnlo/1d22,color=!red,thick=2,lines=2,ps=10

hline,0,col=!gray

legend,[textoidl('\Phi')+'!d<B>,HI!n',textoidl('\Phi')+'!d<B>,LOW!n',textoidl('\Phi')+'!dNET,MF!n'],lines=[0,2,0],color=[!black,!black,!black],thick=[5,2,2],/top,/right


utplot,xarr-mintim,mffluxsigns/1d22,mintim,ps=10,ytit='South '+textoidl('\Phi_{NET}')+' ['+textoidl('\times10^{22}')+' Mx]',pos=[.12,.1,.95,.525],/noerase,yran=[-4,3.99],/ysty,xran=xran,/xsty
oplot,xarr-mintim,mffluxsigns/1d22,color=!blue,thick=2,ps=10

oplot,bsignxarr-mintim,bsignfluxshi/1d22,color=!blue,ps=10,thick=10

oplot,bsignxarr-mintim,bsignfluxslo/1d22,color=!blue,thick=2,lines=2,ps=10
hline,0,col=!gray



closeplotenv
















stop






















end