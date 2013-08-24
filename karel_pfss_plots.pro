;20130819 - make pfss code assimilation plots
pro karel_pfss_plots

root='~/science/projects/Global-Magnetic-Field-Study/'

dpath=root+'./data/'
ppath=root+'./plots/'

restore,dpath+'fluxfile20130805_76.sav',/ver


loadct,0,/sil

tim=(year-1979.)*365.*24.*3600.

setplotenv,file=ppath+'appendix_karel_assymflux.eps',/ps,xs=14,ys=14

!p.multi=[0,1,2]

;plot total, net flux (data in units of 1d18 Mx)

utplot,tim-min(tim),flux/1d4,min(tim),ytit='V2    '+textoidl('\Phi_{TOT} [\times10^{22}')+' Mx]',xtit='',xtickname=strarr(10)+' ',ymarg=[1.5,2],/xsty
utplot,tim-min(tim),netflux/1d4,min(tim),ytit='V2    '+textoidl('\Phi_{NET} [\times10^{22}')+' Mx]',ymarg=[5,-1.5],/xsty

closeplotenv

stop



setplotenv,file=ppath+'appendix_karel_polar.eps',/ps,xs=14,ys=10

!p.multi=0 ;[0,1,2]

;plot polar field strengths

utplot,tim-min(tim),polarflux/1d4,min(tim),ytit='Polar '+textoidl('\Phi_{TOT} [\times10^{22}')+' Mx]',yran=[-5,5],/ysty,/xsty

setcolors,/sys,/sil

hline,0,color=!gray

oplot,tim-min(tim),NORTHPOLARFLUX/1d4,color=!red

oplot,tim-min(tim),SOUTHPOLARFLUX/1d4,color=!blue

legend,['N+S '+textoidl('\Phi_{TOT}'),'N '+textoidl('\Phi_{NET}'),'S '+textoidl('\Phi_{NET}')],lines=[0,0,0],color=[!black,!red,!blue],/bottom,/right

closeplotenv

stop













end