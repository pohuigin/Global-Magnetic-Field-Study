pro smart_avg_prop

restore,'~/science/papers/active_regions_2_cycle/images/smart_sav_1pday.sav',/verb
restore,'~/science/data/cycle23_sav/sunspot_monthly_num_sdic.sav',/ver
;tim_all=anytim(ARSTRUCT_ARR.time)
restore,'~/science/papers/active_regions_2_cycle/images/artim_all.sav',/ver
datapath='~/science/data/cycle23_sav/smart_paper_2/'
paper2path='~/science/papers/active_regions_2_cycle/images/'
ARSTRUCT_ARR=ARSTRUCT_ARR[where(ARSTRUCT_ARR.hglon ge (-60) and ARSTRUCT_ARR.hglon le 60)]
tim_all=anytim(ARSTRUCT_ARR.time)
;tim_all[where(ARSTRUCT_ARR.hglon ge (-60) and ARSTRUCT_ARR.hglon le 60)]

;stop

xbin=.0383562 ;2 week bin = 14./365.
ybin=1.

;Calculate flux imbalance:
fluximb_all=abs((ARSTRUCT_ARR.bfluxpos-ARSTRUCT_ARR.bfluxneg)/ARSTRUCT_ARR.bflux)
fluximbsign_all=(ARSTRUCT_ARR.bfluxpos-ARSTRUCT_ARR.bfluxneg)/ARSTRUCT_ARR.bflux

wsolmin=where((tim_all-min(tim_all)) le 3600.*24.*365.*3.)
wsolmax=where((tim_all-min(tim_all)) gt 3600.*24.*365.*3. and (tim_all-min(tim_all)) le 3600.*24.*365*7)
wsolmin2=where((tim_all-min(tim_all)) gt 3600.*24.*365.*7. and (tim_all-min(tim_all)) le 3600.*24.*365*12)

print,'ALL',anytim(minmax(tim_all),/vms)
print,'MIN',anytim(minmax(tim_all[wsolmin]),/vms)
print,'MAX',anytim(minmax(tim_all[wsolmax]),/vms)
print,'MIN2',anytim(minmax(tim_all[wsolmin2]),/vms)

;NORTH
wnorth=where(ARSTRUCT_ARR.hglat ge 0)
arstructn=ARSTRUCT_ARR[wnorth]

smart_bin_time, arstructn, tim_all[wnorth], n_datastr, n_time, /year, /mean,/smart ;,bin=bin, 
smart_bin_time, arstructn, tim_all[wnorth], n_datastrtot, n_time, /year, /total,/smart

;SOUTH
wsouth=where(ARSTRUCT_ARR.hglat lt 0)
arstructs=ARSTRUCT_ARR[wsouth]
smart_bin_time, arstructs, tim_all[wsouth], s_datastr, s_time, /year, /mean,/smart ;,bin=bin, 
smart_bin_time, arstructs, tim_all[wsouth], s_datastrtot, s_time, /year, /total,/smart

;CYCLE
;smart_bin_time, ARSTRUCT_ARR[wsolmin], tim_all[wsolmin], min_datastr, min_time, /year, /mean,/smart ;,bin=bin, 
;smart_bin_time, ARSTRUCT_ARR[wsolmax], tim_all[wsolmax], max_datastr, max_time, /year, /mean,/smart
;smart_bin_time, ARSTRUCT_ARR[wsolmin2], tim_all[wsolmin2], min2_datastr, min2_time, /year, /mean,/smart

;stop

goto,skipall

;PLOT MEAN PROPERTIES OF ALL REGIONS---------------------------------------------------->
!p.multi=[0,1,3]
setplotenv,/ps,xs=12,ys=12,file=paper2path+'/averageprop/area_flux_all.eps'
setcolors,/sys

utplot,s_time*3600.*24.*365.+anytim(min(tim_all[wsouth]))-anytim(min(tim_all[wnorth])),s_datastr.nar,anytim(min(tim_all[wnorth]),/vms),ytit='# Detections',ps=10,/xsty,chars=4,xtit='',xtickname=strarr(10)+' ',ymargin=[-1,3]
oplot,n_time*3600.*24.*365.,n_datastr.nar,color=!red,ps=10
oplot,sstim-anytim(min(tim_all[wnorth])),ssnum*4d2,color=gray,lines=2

utplot,s_time*3600.*24.*365.+anytim(min(tim_all[wsouth]))-anytim(min(tim_all[wnorth])),s_datastr.area,anytim(min(tim_all[wnorth]),/vms),ytit='AREA [Mm^2]',ps=10,/xsty,chars=4,xtit='',xtickname=strarr(10)+' ',ymargin=[2,1]
oplot,n_time*3600.*24.*365.,n_datastr.area,color=!red,ps=10

utplot,s_time*3600.*24.*365.+anytim(min(tim_all[wsouth]))-anytim(min(tim_all[wnorth])),s_datastr.bflux*1d16,anytim(min(tim_all[wnorth]),/vms),ytit='FLUX [Mx]',ps=10,/xsty,chars=4,ymargin=[4,-2]
oplot,n_time*3600.*24.*365.,n_datastr.bflux*1d16,color=!red,ps=10
closeplotenv & spawn,'convert '+paper2path+'/averageprop/area_flux_all.eps '+paper2path+'/averageprop/area_flux_all.pdf'

skipall:
;PLOT MEAN PROPERTIES OF ACTIVE/MULTIPOLAR REGIONS-------------------------------------->
wnorth=where(ARSTRUCT_ARR.hglat gt 0 and fluximb_all le 0.5) ;take regions with imb. LE half total flux
arstructn=ARSTRUCT_ARR[wnorth]
smart_bin_time, arstructn, tim_all[wnorth], n_datastr, n_time, /year, /mean,/smart ;,bin=bin, 
smart_bin_time, arstructn, tim_all[wnorth], n_datastrtot, n_time, /year, /total,/smart
wsouth=where(ARSTRUCT_ARR.hglat lt 0 and fluximb_all le 0.5)
arstructs=ARSTRUCT_ARR[wsouth]
smart_bin_time, arstructs, tim_all[wsouth], s_datastr, s_time, /year, /mean,/smart ;,bin=bin, 
smart_bin_time, arstructs, tim_all[wsouth], s_datastrtot, s_time, /year, /total,/smart

!p.multi=[0,1,3]
setplotenv,/ps,xs=12,ys=12,file=paper2path+'/averageprop/area_flux_ars.eps'
setcolors,/sys

utplot,s_time*3600.*24.*365.+anytim(min(tim_all[wsouth]))-anytim(min(tim_all[wnorth])),s_datastr.nar,anytim(min(tim_all[wnorth]),/vms),ytit='# Detections',ps=10,/xsty,chars=4,xtit='',xtickname=strarr(10)+' ',ymargin=[-1,3]
oplot,n_time*3600.*24.*365.,n_datastr.nar,color=!red,ps=10
oplot,sstim-anytim(min(tim_all[wnorth])),ssnum*1.3d2,color=gray,lines=2

utplot,s_time*3600.*24.*365.+anytim(min(tim_all[wsouth]))-anytim(min(tim_all[wnorth])),s_datastr.area,anytim(min(tim_all[wnorth]),/vms),ytit='AREA [Mm^2]',ps=10,/xsty,chars=4,xtit='',xtickname=strarr(10)+' ',ymargin=[2,1]
oplot,n_time*3600.*24.*365.,n_datastr.area,color=!red,ps=10

utplot,s_time*3600.*24.*365.+anytim(min(tim_all[wsouth]))-anytim(min(tim_all[wnorth])),s_datastr.bflux*1d16,anytim(min(tim_all[wnorth]),/vms),ytit='FLUX [Mx]',ps=10,/xsty,chars=4,ymargin=[4,-2]
oplot,n_time*3600.*24.*365.,n_datastr.bflux*1d16,color=!red,ps=10
closeplotenv & spawn,'convert '+paper2path+'/averageprop/area_flux_ars.eps '+paper2path+'/averageprop/area_flux_ars.pdf'

;PLOT DISTRIBUTION OF ALL MULTIPOLAR REGIONS-------------------------------------->
wsolmin=where((tim_all-min(tim_all)) le 3600.*24.*365.*3.); and fluximb_all le 0.5)
wsolmax=where((tim_all-min(tim_all)) gt 3600.*24.*365.*3. and (tim_all-min(tim_all)) le 3600.*24.*365*6); and fluximb_all le 0.5)
wsolmin2=where((tim_all-min(tim_all)) gt 3600.*24.*365.*6. and (tim_all-min(tim_all)) le 3600.*24.*365*10); and fluximb_all le 0.5)
!p.multi=[0,2,1]
setplotenv,/ps,xs=20,ys=8,file=paper2path+'/averageprop/area_flux_dist_all.eps'
setcolors,/sys
plot_hist,alog10(ARSTRUCT_ARR.area),xtit='LOG(AREA) [Mm^2]',bin=.05,chars=2,/xsty,/log,xran=[3.5,5.5] ;[where(fluximb_all le 0.5)]
plot_hist,alog10(ARSTRUCT_ARR[wsolmin].area),/oplot,color=!blue,bin=.05,/log
plot_hist,alog10(ARSTRUCT_ARR[wsolmax].area),/oplot,color=!red,bin=.05,/log
plot_hist,alog10(ARSTRUCT_ARR[wsolmin2].area),/oplot,color=!blue,lines=2,bin=.05,/log

!x.tickname=strarr(10)+' ' & !y.tickname=strarr(10)+' '
plot_hist,alog10(ARSTRUCT_ARR[wsolmin].area),xtit=' ',ytit=' ',color=!blue,bin=.05,/log,/xsty,/ysty,/noerase,chars=2,xran=[3.5,5.5];,xtickname=strarr(10)+' ',ytickname=' '
plot_hist,alog10(ARSTRUCT_ARR[wsolmax].area),xtit=' ',ytit=' ',color=!red,bin=.05,/log,/xsty,/ysty,/noerase,chars=2,xran=[3.5,5.5];,xtickname=strarr(10)+' ',ytickname=' '
plot_hist,alog10(ARSTRUCT_ARR[wsolmin2].area),xtit=' ',ytit=' ',color=!blue,lines=2,bin=.05,/log,/xsty,/ysty,/noerase,chars=2,xran=[3.5,5.5];,xtickname=strarr(10)+' ',ytickname=' '
!x.tickname='' & !y.tickname=''
plot_hist,alog10(ARSTRUCT_ARR.area),ytit=' ',xtit='LOG(AREA) [Mm^2]',bin=.05,chars=2,/xsty,/ysty,/log,/noerase,xran=[3.5,5.5];[where(fluximb_all le 0.5)]
closeplotenv & spawn,'convert '+paper2path+'/averageprop/area_flux_dist_all.eps '+paper2path+'/averageprop/area_flux_dist_all.pdf'

!p.multi=[0,2,1]
setplotenv,/ps,xs=20,ys=8,file=paper2path+'/averageprop/flux_dist_all.eps'
setcolors,/sys
plot_hist,alog10(ARSTRUCT_ARR.bflux*1d16),xtit='LOG(FLUX) [Mx]',bin=.05,chars=2,/xsty,/log,xran=[22,24],yran=[1d,1d5],/ysty;[fluximb_all le 0.5]
plot_hist,alog10(ARSTRUCT_ARR[wsolmin].bflux*1d16),/oplot,color=!blue,bin=.05,/log,xran=[22,24],yran=[1d,1d5],/ysty
plot_hist,alog10(ARSTRUCT_ARR[wsolmax].bflux*1d16),/oplot,color=!red,bin=.05,/log,xran=[22,24],yran=[1d,1d5],/ysty
plot_hist,alog10(ARSTRUCT_ARR[wsolmin2].bflux*1d16),/oplot,color=!blue,lines=2,bin=.05,/log,xran=[22,24],yran=[1d,1d5],/ysty

!x.tickname=strarr(10)+' ' & !y.tickname=strarr(10)+' '
plot_hist,alog10(ARSTRUCT_ARR[wsolmin].bflux*1d16),xtit=' ',ytit=' ',color=!blue,bin=.05,/log,/xsty,/ysty,/noerase,chars=2,xran=[22,24];,xtickname=strarr(10)+' ',ytickname=' '
plot_hist,alog10(ARSTRUCT_ARR[wsolmax].bflux*1d16),xtit=' ',ytit=' ',color=!red,bin=.05,/log,/xsty,/ysty,/noerase,chars=2,xran=[22,24];,xtickname=strarr(10)+' ',ytickname=' '
plot_hist,alog10(ARSTRUCT_ARR[wsolmin2].bflux*1d16),xtit=' ',ytit=' ',color=!blue,lines=2,bin=.05,/log,/xsty,/ysty,/noerase,chars=2,xran=[22,24];,xtickname=strarr(10)+' ',ytickname=' '
!x.tickname='' & !y.tickname=''
plot_hist,alog10(ARSTRUCT_ARR.bflux*1d16),ytit='Normalized to Max',xtit='LOG(FLUX) [Mx]',bin=.05,chars=2,/xsty,/ysty,/log,/noerase,xran=[22,24];[fluximb_all le 0.5]
closeplotenv & spawn,'convert '+paper2path+'/averageprop/flux_dist_all.eps '+paper2path+'/averageprop/flux_dist_all.pdf'

;!p.multi=[0,2,1]
;setplotenv,/ps,xs=20,ys=8,file=paper2path+'/averageprop/rval_dist_all.eps'
;setcolors,/sys
;plot_hist,alog10(ARSTRUCT_ARR.nlstr.rval),xtit='LOG(FLUX) [Mx]',bin=.1,chars=2,/xsty,/log
;plot_hist,alog10(ARSTRUCT_ARR[wsolmin].nlstr.rval),/oplot,color=!blue,bin=.1,/log
;plot_hist,alog10(ARSTRUCT_ARR[wsolmax].nlstr.rval),/oplot,color=!red,bin=.1,/log
;plot_hist,alog10(ARSTRUCT_ARR[wsolmin2].nlstr.rval),/oplot,color=!blue,lines=2,bin=.1,/log

;!x.tickname=strarr(10)+' ' & !y.tickname=strarr(10)+' '
;plot_hist,alog10(ARSTRUCT_ARR[wsolmin].nlstr.rval),xtit=' ',ytit=' ',color=!blue,bin=.1,/log,/xsty,/ysty,/noerase,chars=2,xran=[10,16];,xtickname=strarr(10)+' ',ytickname=' '
;plot_hist,alog10(ARSTRUCT_ARR[wsolmax].nlstr.rval),xtit=' ',ytit=' ',color=!red,bin=.1,/log,/xsty,/ysty,/noerase,chars=2,xran=[10,16];,xtickname=strarr(10)+' ',ytickname=' '
;plot_hist,alog10(ARSTRUCT_ARR[wsolmin2].nlstr.rval),xtit=' ',ytit=' ',color=!blue,lines=2,bin=.1,/log,/xsty,/ysty,/noerase,chars=2,xran=[10,16];,xtickname=strarr(10)+' ',ytickname=' '
;!x.tickname='' & !y.tickname=''
;plot_hist,alog10(ARSTRUCT_ARR.nlstr.rval),ytit='Normalized to Max',xtit='LOG(R) [Mx]',bin=.1,chars=2,/xsty,/ysty,/log,/noerase,xran=[10,16]
;closeplotenv & spawn,'convert '+paper2path+'/averageprop/rval_dist_all.eps '+paper2path+'/averageprop/rval_dist_all.pdf'

stop






;PLOT DISTRIBUTIONS FOR MIN MAX CYCLE---------------------------------->
yharea=histogram(alog10(ARSTRUCT_ARR.area),locat=xharea)
bm=linfit(xharea[where(xharea ge 3.3)],yharea[where(xharea ge 3.3)])
oplot,xharea[where(xharea ge 3.3)],xharea[where(xharea ge 3.3)]*bm[1]+bm[0],color=!red

yharea=histogram(ARSTRUCT_ARR.area,locat=xharea,bin=1d3)
param=mpfitexpr('P[0]*(P[1])^(X)', xharea[where(xharea ge 2d3 and xharea le 7d4)], yharea[where(xharea ge 2d3 and xharea le 7d4)],sy,[1d,1d])

stop

;FIT POWERLAW (diff for N and south?)




;PLOT DISTRIBUTIONS FOR DIFFERENT TIME PERIODS

!p.multi=[0,1,3]
setplotenv,/ps,xs=12,ys=12,file=paper2path+'/averageprop/area_flux_dist_all.eps'
setcolors,/sys





;PLOT POWER LAW ALPHA AS A FUNCTION OF YEAR IN CYCLE

stop

n_time,n_datastr.flux
n_time,n_datastr.netflux
n_time,n_datastr.fluximb






end