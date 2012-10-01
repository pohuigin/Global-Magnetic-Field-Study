pro smart_plot_ar_flr_carrington

restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday.sav',/ver
restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_tim.sav',/ver ;checked! in units of years
restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_carrarr.sav',/ver

restore,'~/science/data/goes_event_list/goes_flare_struct.sav',/ver

;check carringlon longitude calculation
;carlonarr2=tim2carr(arstr_arr.time)+arstr_arr.hglon
;wclt0=where(carlonarr2 lt 0)
;if wclt0[0] ne -1 then carlonarr2[wclt0]=carlonarr2[wclt0]+360.
;wcge360=where(carlonarr2 ge 360)
;if wcge360[0] ne -1 then carlonarr2[wcge360]=carlonarr2[wcge360]-360.
;stop

datapath='~/science/data/cycle23_sav/smart_paper_2/'

;Take flares between ar times where theres data
ARSTR_ARR=ARSTR_ARR[where(ARSTR_ARR.hglat ne 10000)]
TIM1PD=TIM1PD[where(ARSTR_ARR.hglat ne 10000)]
CARLONARR=CARLONARR[where(ARSTR_ARR.hglat ne 10000)]
wars=where((FLRPOSSTR.tim/3600./24./365.-mintim) ge min(tim1pd) and (FLRPOSSTR.tim/3600./24./365.-mintim) le max(tim1pd))
astrpos=astrpos[wars]
FLRPOSSTR=FLRPOSSTR[wars]

setplotenv,/xwin & !p.color=0 & !p.background=255
window,8,xs=1000,ys=700

;PLOT Carrington Longitude flares complexity----------------------------------------->
plotsym,0,.5,/fill
setcolors,/sys
!p.multi=0

plot,tim1pd,carlonarr,ps=8,yran=[0,360],xran=[0,11.2],/xsty,/ysty,xtit='years since 19970101',ytit='carrington longitude',tit='R=10^22mx G=15k WLSG B=MX Flares'
oplot,tim1pd,carlonarr,ps=8,color=150

plotsym,0,.7,/fill
wgt22=where(arstr_arr.bflux*1d16 gt 1d22)
oplot,tim1pd[wgt22],carlonarr[wgt22],ps=8,color=!red

wwlsg=where(arstr_arr.nlstr.wlsg gt 15000)
oplot,tim1pd[wwlsg],carlonarr[wwlsg],ps=8,color=!green
wle60=where(arstr_arr[wwlsg].hglon le 60.)
oplot,tim1pd[wwlsg[wle60]],carlonarr[wwlsg[wle60]],ps=8,color=!forest

plotsym,0,.7
wmx=where(strmid(ASTRPOS.GOES_CLASS,0,1) eq 'M' or strmid(ASTRPOS.GOES_CLASS,0,1) eq 'X')
window_capture,file=datapath+'carlon_wlsg_mxflares_plot0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,FLRPOSSTR[wmx].carlon,ps=8,color=0
window_capture,file=datapath+'carlon_wlsg_mxflares_plot1',/png

stop

;PLOT HOWARD diff_rot Longitude flares complexity----------------------------------------->
plotsym,0,.5,/fill
setcolors,/sys
!p.multi=0

howardlon=(diff_rot(mintim/2.-365.*tim1pd,ARSTR_ARR.hglat) mod 360.)+ARSTR_ARR.hglon
howardlon[where(howardlon lt 0)]=howardlon[where(howardlon lt 0)]+360.
flrhowardlon=abs((diff_rot(-365.*(FLRPOSSTR.tim/3600./24./365.-mintim),FLRPOSSTR.hgpos[1,*]))+FLRPOSSTR.hgpos[0,*]) mod 360.)

plot,tim1pd,howardlon,ps=8,yran=[0,360],xran=[0,11.2],/xsty,/ysty,xtit='years since 19970101',ytit='carrington longitude',tit='R=10^22mx G=15k WLSG B=MX Flares'
oplot,tim1pd,howardlon,ps=8,color=150

plotsym,0,.7,/fill
oplot,tim1pd[wgt22],howardlon[wgt22],ps=8,color=!red
oplot,tim1pd[wwlsg],howardlon[wwlsg],ps=8,color=!green
oplot,tim1pd[wwlsg[wle60]],howardlon[wwlsg[wle60]],ps=8,color=!forest
plotsym,0,.7
window_capture,file=datapath+'carlon_howard_wlsg_mxflares_plot0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,flrhowardlon[wmx],ps=8,color=0
window_capture,file=datapath+'carlon_howard_wlsg_mxflares_plot1',/png

stop

;PLOT ALLEN diff_rot Longitude flares complexity----------------------------------------->









end