pro smart_polar_plots

datapath='~/science/data/cycle23_sav/smart_paper_2/'

plotpath='~/science/talks/sgm_20100618/plots/'

restore,'~/science/data/smart_sf_smart_compare/flare_data_sf.sav',/ver

restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday.sav',/ver

restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_tim.sav',/ver

restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_carrarr.sav',/ver

;restore polar data
read_dipole_data, polestruct, /res

restore,'~/science/data/cycle23_sav/sunspot_monthly_num_sdic.sav',/ver

restore,'~/science/talks/arg_retreat_20100514/smart_cycle.sav',/ver

restore,'~/science/data/goes_event_list/goes_flare_struct.sav',/ver

xmarginorig=!x.margin

setplotenv,/xwin
!p.multi=0
!p.color=0
!p.background=255

window,8,xs=800,ys=500

;-------------------------------------------------------->

plot,SSTIM/3600./24./365.+1979.,ssnum,ytit='SIDC Sunspot #',thick=3,xtit='Years'
window_capture,file=plotpath+'sidc_sunspot_num',/png

;-------------------------------------------------------->

!x.margin=[7,7]
setcolors,/sys,/quiet
plot,polestruct.tim/3600./24./365.+1979.,polestruct.filt_nth*10.,color=!red, xr=[1996,2010],/xsty,yr=[-1100,1100],/ysty,xticklen=.0001, yticklen=.0001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' '
oplot,polestruct.tim/3600./24./365.+1979.,polestruct.filt_sth*10.,color=!blue
hline,0,lines=2
plot,SSTIM/3600./24./365.+1979.,ssnum,ytit='SIDC Sunspot #',thick=3,xtit='Years', xr=[1996,2010],/xsty,/noerase,yr=[0,140],/ysty
axis,yaxis=1,color=!white,yran=[0,140],/ysty
axis,yaxis=1,color=!red,yran=[-1100,1100],/ysty,ytit='B Poles [mG]'
xyouts,.73,.8,'North',color=!red,/norm
xyouts,.80,.8,'South',color=!blue,/norm
window_capture,file=plotpath+'sidc_vs_poles',/png

;-------------------------------------------------------->

wa=where(strmid(smflrarr.class,0,1) eq 'A')
wb=where(strmid(smflrarr.class,0,1) eq 'B')
wc=where(strmid(smflrarr.class,0,1) eq 'C')
wm=where(strmid(smflrarr.class,0,1) eq 'M')
wx=where(strmid(smflrarr.class,0,1) eq 'X')

wgt23=where(arstr_arr.bflux*1d16 gt 1d23)
wgt22=where(arstr_arr.bflux*1d16 gt 1d22)
wgt21=where(arstr_arr.bflux*1d16 ge 1d21)
wlt21=where(arstr_arr.bflux*1d16 lt 1d21)

fluximbsign=(arstr_arr.bfluxpos-arstr_arr.bfluxneg)/arstr_arr.bflux
w0tp3=where(fluximbsign ge 0 and fluximbsign le .3)
wp3tp6=where(fluximbsign gt .3 and fluximbsign le .6)
wp6tp1=where(fluximbsign gt .6 and fluximbsign le 1)
w0tn3=where(fluximbsign lt 0 and fluximbsign ge -.3)
wn3tn6=where(fluximbsign lt -.3 and fluximbsign lt -.6)
wn6tn1=where(fluximbsign lt -.6 and fluximbsign ge -1)

wp8tp1=where(fluximbsign gt .8 and fluximbsign le 1)
wn8tn1=where(fluximbsign lt -.8 and fluximbsign ge -1)
wn8tp8=where(fluximbsign ge -.8 and fluximbsign le .8)

;-------------------------------------------------------->
;plor db_poles/dt with smart AR #

window,8,xs=700,ys=700
!p.multi=[0,1,2]
plotsym,0,.6,/fill

smhisty=histogram(tim1pd,loc=smhistx,bin=.1)
;plot_hist,tim1pd,ytit='SMART AR #',bin=.1,xran=[0,11.2],/xsty,charsize=1.8,xtit='Years since Feb 3, 1997'
plot,smhistx,smhisty/36.5,xran=[0,11.2],/xsty,ps=10,charsize=1.8, ytit='SMART AR #',pos=[.1,.6,.9,.9];,xtit='Years since Feb 3, 1997'
totfluxarrsmth=smooth(totfluxarr,36.5*2.)
totfluxnorthsmth=smooth(totfluxnorth,36.5*2.)
totfluxsouthsmth=smooth(totfluxsouth,36.5*2.)
oplot,utim,totfluxarrsmth/max(totfluxarrsmth)*max(smhisty)/36.5,color=150,ps=8
oplot,utim,totfluxnorthsmth/max(totfluxarrsmth)*max(smhisty)/36.5,color=!red,ps=8
oplot,utim,totfluxsouthsmth/max(totfluxarrsmth)*max(smhisty)/36.5,color=!blue,ps=8

;plot,(polestruct.tim-anytim('2-feb-1997'))/3600./24./365.,polestruct.filt_nth*10.,color=!red,/xsty,yr=[-1100,1100],xr=[0,11.2],/ysty
;oplot,(polestruct.tim-anytim('2-feb-1997'))/3600./24./365.,polestruct.filt_sth*10.,color=!blue
;hline,0,lines=2

plot,(polestruct.tim-anytim('2-feb-1997'))/3600./24./365.,deriv(polestruct.tim,polestruct.filt_nth),color=!red,/xsty,xr=[0,11.2]
oplot,(polestruct.tim-anytim('2-feb-1997'))/3600./24./365.,deriv(polestruct.tim,polestruct.filt_sth),color=!blue
hline,0,lines=2

;mintim

;-------------------------------------------------------->
;Plot hemispheric AR # vs hemispheric pole strength

!p.multi=[0,1,2]

npoleorigx=((polestruct.tim-anytim('2-feb-1997'))/3600./24./365.)
npoleorigy=deriv(polestruct.tim,polestruct.filt_nth)
wgoodnpole=where(npoleorigx ge min(utim))
npolex=npoleorigx[wgoodnpole]
npoley=npoleorigy[wgoodnpole]

wgoodnsmar=where(utim le max(npolex))
nsmarx=utim[wgoodnsmar]
nsmary=totfluxnorthsmth[wgoodnsmar]
npoleyint=interpol(npoley,npolex,nsmarx,/quad)

plot,npolex,npoley,tit='North Hemisphere AR Flux and Pole Field'
hline,0
oplot,nsmarx,npoleyint,color=!red
oplot,nsmarx,nsmary*1d-29,color=!blue
;plot,nsmary,npoleyint
;oplot,nsmary,npoleyint,ps=4,color=!red

spoleorigx=((polestruct.tim-anytim('2-feb-1997'))/3600./24./365.)
spoleorigy=deriv(polestruct.tim,polestruct.filt_sth)
wgoodspole=where(spoleorigx ge min(utim))
spolex=spoleorigx[wgoodspole]
spoley=spoleorigy[wgoodspole]

wgoodssmar=where(utim le max(spolex))
ssmarx=utim[wgoodssmar]
ssmary=totfluxsouthsmth[wgoodssmar]
spoleyint=interpol(spoley,spolex,ssmarx,/quad)

plot,spolex,spoley,tit='South Hemisphere AR Flux and Pole Field'
hline,0
oplot,ssmarx,spoleyint,color=!red
oplot,ssmarx,ssmary*1d-29,color=!blue

stop

;-------------------------------------------------------->
;Plot butterfly diagram

window,8,xs=1000,ys=700
!p.multi=0

plot,tim1pd,arstr_arr.hglat,yran=[-60,60],ps=3,/xsty,xran=[0,11.2],xtit='Years since '+anytim(arstr_arr[0].time,/vms,/date),ytit='Latitude',pos=[.1,.15,.95,.9]
oplot,tim1pd[wlt21],arstr_arr[wlt21].hglat,ps=8,color=150
oplot,tim1pd[wgt21],arstr_arr[wgt21].hglat,ps=8,color=100
oplot,tim1pd[wgt22],arstr_arr[wgt22].hglat,ps=8,color=50
oplot,tim1pd[wgt23],arstr_arr[wgt23].hglat,ps=8,color=0
window_capture,file=datapath+'hglat_butterfly_bflux',/png

plotsym,0,.5
wbc=where(strmid(ASTRPOS.GOES_CLASS,0,1) eq 'B' or strmid(ASTRPOS.GOES_CLASS,0,1) eq 'C')
oplot,FLRPOSSTR[wbc].tim/(3600.*24.*365)-mintim,(FLRPOSSTR[wbc].hgpos)[1,*],ps=8,color=!green
wmx=where(strmid(ASTRPOS.GOES_CLASS,0,1) eq 'M' or strmid(ASTRPOS.GOES_CLASS,0,1) eq 'X')
oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,(FLRPOSSTR[wmx].hgpos)[1,*],ps=8,color=!red
window_capture,file=datapath+'hglat_butterfly_bflux_flares',/png
stop

;-------------------------------------------------------->
;Plot butterfly carrington longitude

window,8,xs=1000,ys=600
!p.multi=0

;carlonarr=arstr_arr.carlon
;carlonarr[where(carlonarr lt 0)]=carlonarr[where(carlonarr lt 0)]+360.

;Restored carlonarr at beginning
	;carlonarr=tim2carr(arstr_arr.time)
	;carlonarr=carlonarr+arstr_arr.hglon
	;carlonarr[where(carlonarr ge 360)]=carlonarr[where(carlonarr ge 360)]-360.
	;carlonarr[where(carlonarr lt 0)]=carlonarr[where(carlonarr lt 0)]+360.
wnot0=where(arstr_arr.hglat ne 10000)

;plain
;original carrlon
;plot,(tim1pd)[wnot0],(arstr_arr.carlon)[wnot0],/ysty,yran=[0,360],ps=3,/xsty,xran=[0,11.2],xtit='Years since '+anytim(arstr_arr[0].time,/vms,/date),ytit='Longitude',pos=[.1,.15,.95,.9]
;oplot,(tim1pd)[wlt21],(arstr_arr.carlon)[wlt21],ps=8,color=150
;oplot,tim1pd[wgt21],(arstr_arr.carlon)[wgt21],ps=8,color=100
;oplot,tim1pd[wgt22],(arstr_arr.carlon)[wgt22],ps=8,color=50
;oplot,tim1pd[wgt23],(arstr_arr.carlon)[wgt23],ps=8,color=0
;stop
plot,tim1pd[wnot0],carlonarr[wnot0],/ysty,yran=[0,360],ps=3,/xsty,xran=[0,11.2],xtit='Years since '+anytim(arstr_arr[0].time,/vms,/date),ytit='Longitude',pos=[.1,.15,.95,.9]
oplot,tim1pd[wlt21],carlonarr[wlt21],ps=8,color=150
oplot,tim1pd[wgt21],carlonarr[wgt21],ps=8,color=100
oplot,tim1pd[wgt22],carlonarr[wgt22],ps=8,color=50
oplot,tim1pd[wgt23],carlonarr[wgt23],ps=8,color=0

wnorth=where(arstr_arr.hglat gt 0)
wnorthgt22=wnorth[where((arstr_arr.bflux)[wnorth]*1d16 gt 1d22)]
wsouth=where(arstr_arr.hglon le 0)
wsouthgt22=wsouth[where((arstr_arr.bflux)[wsouth]*1d16 gt 1d22)]
oplot,tim1pd[wnorthgt22],carlonarr[wnorthgt22],ps=8,color=!red
oplot,tim1pd[wsouthgt22],carlonarr[wsouthgt22],ps=8,color=!blue

stop

;zoomed
plot,tim1pd[wnot0],carlonarr[wnot0],/ysty,yran=[0,360],ps=3,/xsty,xran=[2.5,4],xtit='Years since '+anytim(arstr_arr[0].time,/vms,/date),ytit='Longitude',pos=[.1,.15,.95,.9]
oplot,tim1pd[wlt21],carlonarr[wlt21],ps=8,color=150
oplot,tim1pd[wgt21],carlonarr[wgt21],ps=8,color=100
oplot,tim1pd[wgt22],carlonarr[wgt22],ps=8,color=50
oplot,tim1pd[wgt23],carlonarr[wgt23],ps=8,color=0
oplot,tim1pd[wnorthgt22],carlonarr[wnorthgt22],ps=8,color=!red
oplot,tim1pd[wsouthgt22],carlonarr[wsouthgt22],ps=8,color=!blue

stop

;-------------------------------------------------------->

;over plot large regions integrated over lat,lon on all regions over time. see if there is a phase shift. all large regions occur at beginning of burst of emergence

;plot dphi/dt over time, compare to polar db/dt or polar strengths in general 

;compare global flux in unipolar+ and unipolar- in north and south hemispheres.
;compare to polar fields

;compare number of large regions on disk to dipolarity ratio, also unipolar regions, solar cycle phase, etc.





!x.margin=xmarginorig

stop
end