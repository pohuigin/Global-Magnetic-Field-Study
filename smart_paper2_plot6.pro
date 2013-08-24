;@smart_paper2_plot6.pro
;Make plot 6 and 6_5 for the paper

;routine to test plot the v1 and v2 PFSS runs and compare the a,b, and C's...

datapath='./data/'
plotpath='./plots/20110330/'

restore,datapath+'smart_cycle_edge_linfits_20110330.sav' ;to get WNCENT and WSCENT
restore,datapath+'cycle23_butterfly_flux_maps_final.sav',/ver
restore,datapath+'butterfly_diags_lat_27day_bin.sav',/ver
;restore,datapath+'pfss_monthly_phiab_plot.sav',/ver
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

.r sphere_harmonic

;VERSION 2-------------->

harmfile=datapath+'pfss_monthly_phiab_coeff_v2.sav'
restore,harmfile,/ver
A_lm=PHIATARR
B_lm=PHIBTARR
abtimarr=anytim(date)
ntime=n_elements(date)

v2coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1, /use_phi,outa=a00, outb=b00)
v2coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1, /use_phi,outa=a01, outb=b01)
v2coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1, /use_phi,outa=a02, outb=b02)
v2coeff03=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=1, /use_phi,outa=a03, outb=b03)
v2coeff04=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=1, /use_phi,outa=a04, outb=b04)
v2coeff05=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=1, /use_phi,outa=a05, outb=b05)

v2coeff00r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1.75, /use_phi,outa=a00, outb=b00)
v2coeff01r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1.75, /use_phi,outa=a01, outb=b01)
v2coeff02r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1.75, /use_phi,outa=a02, outb=b02)
v2coeff03r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=1.75, /use_phi,outa=a03, outb=b03)
v2coeff04r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=1.75, /use_phi,outa=a04, outb=b04)
v2coeff05r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=1.75, /use_phi,outa=a05, outb=b05)

;yaxis value to over plot and indicate the COEFF of greatest magnitude
plevel=1.8

;Smooth by ~3 rotations and determine which is strongest at each time
v2coeff02sm=abs(smooth(v2coeff02,5))
v2coeff04sm=abs(smooth(v2coeff04,5))
w02=where(v2coeff02sm gt v2coeff04sm)
w04=where(v2coeff04sm gt v2coeff02sm)

c02gtarrv2=fltarr(ntime)*1./0. & c02gtarrv2[w02]=plevel
c04gtarrv2=fltarr(ntime)*1./0. & c04gtarrv2[w04]=plevel

v2coeff01sm=abs(smooth(v2coeff01,5))
v2coeff03sm=abs(smooth(v2coeff03,5))
v2coeff05sm=abs(smooth(v2coeff05,5))
w01=where(v2coeff01sm gt v2coeff03sm and v2coeff01sm gt v2coeff05sm)
w03=where(v2coeff03sm gt v2coeff01sm and v2coeff03sm gt v2coeff05sm)
w05=where(v2coeff05sm gt v2coeff01sm and v2coeff05sm gt v2coeff03sm)

c01gtarrv2=fltarr(ntime)*1./0. & c01gtarrv2[w01]=plevel
c03gtarrv2=fltarr(ntime)*1./0. & c03gtarrv2[w03]=plevel
c05gtarrv2=fltarr(ntime)*1./0. & c05gtarrv2[w05]=plevel



;THEN PLOT IT TO ON PLOT_6 to INDICATE WHICH COEFF IS GTest AT EACH TIME!


;stop


yr=(abtimarr-anytim('01-jan-1996 00:00:00'))/3600./24./365.

;window,xs=800,ys=800
;!p.multi=[0,1,3]
;!y.margin=[3,0]
;!p.charsize=2
;!y.range=[-2,2]
;!x.style=1

;plot,yr,a00,ytit='m0,l0',ps=10
;oplot,yr,b00,lines=2,ps=10                
;plot,yr,a01,ytit='m0,l1',ps=10
;oplot,yr,b01,lines=2,ps=10                
;plot,yr,a02,ytit='m0,l2',ps=10
;oplot,yr,b02,lines=2,ps=10                

;imgpath='~/science/papers/active_regions_2_cycle/images/pfss_test_v1v2/'
;window_capture,file=imgpath+'test_v2_ab'






;PLOT 6. - plot polar fields and dipole sphere harmonics---------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_6_polarb_sphereharm_v2.eps'
setcolors,/sys,/silent,/quiet
;mintim=min(tflux27dseries)
mintim=min(ABTIMARR)

linecol=0

;polar plot
setcolors,/sys
utplot,timrebin-mintim,npole6065,mintim,yran=[-10,10], pos=[.12,.66,.95,.95],ps=10,xran=minmax(timrebin)-mintim,/xsty,ytit='Polar Fields [G]',xtit='',xtickname=strarr(10)+' '
oplot,timrebin-mintim,npole6065,color=!red,ps=10
oplot,timrebin-mintim,npole6570,color=!red,ps=10,lines=2
;oplot,timrebin-mintim,npole7075,color=!red,ps=10
oplot,timrebin-mintim,spole6065,color=!blue,ps=10
oplot,timrebin-mintim,spole6570,color=!blue,ps=10,lines=2
;oplot,timrebin-mintim,spole7075,color=!blue,ps=10
hline,0,color=!gray

legend,['Lat. 60-65','Lat. 65-70'],/top,/right,color=[linecol,linecol],lines=[0,2]


;spherical harmonic plot
;monopole
utplot,ABTIMARR-mintim,reform(v2COEFF02),mintim,pos=[.12,.37,.95,.66],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-2,1.99],ytit='Harmonic Coeff. (R=R!dSun!n)',xtit='',xtickname=strarr(10)+' '
;3 bands total
oplot,ABTIMARR-mintim,reform(v2COEFF00),ps=10,color=!gray
;5 bands total
oplot,ABTIMARR-mintim,reform(v2COEFF04),ps=10,lines=2
hline,0,color=!gray
oplot,ABTIMARR-mintim,c02gtarrv2,thick=10
oplot,ABTIMARR-mintim,c04gtarrv2,thick=5,lines=2


legend,[textoidl('C_{0,0}'),textoidl('C_{2,0}'),textoidl('C_{4,0}')],/bottom,/right,color=[!gray,linecol,linecol],lines=[0,0,2]

;pure dipole
utplot,ABTIMARR-mintim,reform(v2COEFF01),mintim,pos=[.12,.1,.95,.37],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-2,1.99],ytit='Harmonic Coeff. (R=R!dSun!n)'
;two bands in each hemisphere
oplot,ABTIMARR-mintim,reform(v2COEFF03),ps=10,color=!gray
;three bands in each hemisphere (AR, decay, pole)
oplot,ABTIMARR-mintim,reform(v2COEFF05),ps=10,lines=2

;oplot,ABTIMARR-mintim,reform(v2COEFF02),ps=10,color=!red ;non hemisherically symmetric
;oplot,ABTIMARR-mintim,reform(v2COEFF04),ps=10,color=!green ;non hemisherically symmetric
hline,0,color=!gray
oplot,ABTIMARR-mintim,c01gtarrv2,thick=10
oplot,ABTIMARR-mintim,c03gtarrv2,thick=10,color=!gray
oplot,ABTIMARR-mintim,c05gtarrv2,thick=5,lines=2

legend,[textoidl('C_{1,0}'),textoidl('C_{3,0}'),textoidl('C_{5,0}')],/bottom,/right,color=[linecol,!gray,linecol],lines=[0,0,2]

closeplotenv

;----------------------------------------------------------------------------->





;yaxis value to over plot and indicate the COEFF of greatest magnitude
plevelr175=0.45

;Smooth by ~3 rotations and determine which is strongest at each time
v2coeff02r175sm=abs(smooth(v2coeff02r175,5))
v2coeff04r175sm=abs(smooth(v2coeff04r175,5))
w02=where(v2coeff02r175sm gt v2coeff04r175sm)
w04=where(v2coeff04r175sm gt v2coeff02r175sm)

c02gtarrv2r175=fltarr(ntime)*1./0. & c02gtarrv2r175[w02]=plevelr175
c04gtarrv2r175=fltarr(ntime)*1./0. & c04gtarrv2r175[w04]=plevelr175

v2coeff01r175sm=abs(smooth(v2coeff01r175,5))
v2coeff03r175sm=abs(smooth(v2coeff03r175,5))
v2coeff05r175sm=abs(smooth(v2coeff05r175,5))
w01=where(v2coeff01r175sm gt v2coeff03r175sm and v2coeff01r175sm gt v2coeff05r175sm)
w03=where(v2coeff03r175sm gt v2coeff01r175sm and v2coeff03r175sm gt v2coeff05r175sm)
w05=where(v2coeff05r175sm gt v2coeff01r175sm and v2coeff05r175sm gt v2coeff03r175sm)

c01gtarrv2r175=fltarr(ntime)*1./0. & c01gtarrv2r175[w01]=plevelr175
c03gtarrv2r175=fltarr(ntime)*1./0. & c03gtarrv2r175[w03]=plevelr175
c05gtarrv2r175=fltarr(ntime)*1./0. & c05gtarrv2r175[w05]=plevelr175





;PLOT 6.1 - plot sphere harmonics for R=1.75---------------------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_6_1_polarb_sphereharm_r1_75_v2.eps'
setcolors,/sys,/silent,/quiet
;mintim=min(tflux27dseries)
mintim=min(ABTIMARR)

linecol=0

;spherical harmonic plot
;monopole
utplot,ABTIMARR-mintim,reform(v2coeff02r175),mintim,pos=[.12,.55,.95,.95],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-.5,.5],/ysty,ytit='Harmonic Coeff. (R=1.75R!dSun!n)',xtit='',xtickname=strarr(10)+' '
;3 bands total
oplot,ABTIMARR-mintim,reform(v2coeff00r175),ps=10,color=!gray
;5 bands total
oplot,ABTIMARR-mintim,reform(v2coeff04r175),ps=10,lines=2
hline,0,color=!gray
oplot,ABTIMARR-mintim,c02gtarrv2r175,thick=10
oplot,ABTIMARR-mintim,c04gtarrv2r175,thick=5,lines=2

legend,[textoidl('C_{0,0}'),textoidl('C_{2,0}'),textoidl('C_{4,0}')],/bottom,/right,color=[!gray,linecol,linecol],lines=[0,0,2]

;pure dipole
utplot,ABTIMARR-mintim,reform(v2coeff01r175),mintim,pos=[.12,.1,.95,.55],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-.5,.5],/ysty,ytit='Harmonic Coeff. (R=1.75R!dSun!n)'
;two bands in each hemisphere
oplot,ABTIMARR-mintim,reform(v2coeff03r175),ps=10,color=!gray
;three bands in each hemisphere (AR, decay, pole)
oplot,ABTIMARR-mintim,reform(v2coeff05r175),ps=10,lines=2

;oplot,ABTIMARR-mintim,reform(COEFF02),ps=10,color=!red ;non hemisherically symmetric
;oplot,ABTIMARR-mintim,reform(COEFF04),ps=10,color=!green ;non hemisherically symmetric
hline,0,color=!gray
oplot,ABTIMARR-mintim,c01gtarrv2r175,thick=10
oplot,ABTIMARR-mintim,c03gtarrv2r175,thick=10,color=!gray
oplot,ABTIMARR-mintim,c05gtarrv2r175,thick=5,lines=2

legend,[textoidl('C_{1,0}'),textoidl('C_{3,0}'),textoidl('C_{5,0}')],/bottom,/right,color=[linecol,!gray,linecol],lines=[0,0,2]

closeplotenv













;VERSION 1------------------>

harmfile=datapath+'pfss_monthly_phiab_coeff.sav'
restore,harmfile,/ver
A_lm=PHIATARR
B_lm=PHIBTARR
abtimarr=anytim(date)

v1coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1, /use_phi,outa=a00, outb=b00)
v1coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1, /use_phi,outa=a01, outb=b01)
v1coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1, /use_phi,outa=a02, outb=b02)
v1coeff03=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=1, /use_phi,outa=a03, outb=b03)
v1coeff04=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=1, /use_phi,outa=a04, outb=b04)
v1coeff05=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=1, /use_phi,outa=a05, outb=b05)

v1coeff00r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1.75, /use_phi,outa=a00, outb=b00)
v1coeff01r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1.75, /use_phi,outa=a01, outb=b01)
v1coeff02r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1.75, /use_phi,outa=a02, outb=b02)
v1coeff03r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=1.75, /use_phi,outa=a03, outb=b03)
v1coeff04r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=1.75, /use_phi,outa=a04, outb=b04)
v1coeff05r175=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=1.75, /use_phi,outa=a05, outb=b05)

;yaxis value to over plot and indicate the COEFF of greatest magnitude
plevel=1.8

;Smooth by ~3 rotations and determine which is strongest at each time
v1coeff02sm=abs(smooth(v1coeff02,5))
v1coeff04sm=abs(smooth(v1coeff04,5))
w02=where(v1coeff02sm gt v1coeff04sm)
w04=where(v1coeff04sm gt v1coeff02sm)

c02gtarrv1=fltarr(ntime)*1./0. & c02gtarrv1[w02]=plevel
c04gtarrv1=fltarr(ntime)*1./0. & c04gtarrv1[w04]=plevel

v1coeff01sm=abs(smooth(v1coeff01,5))
v1coeff03sm=abs(smooth(v1coeff03,5))
v1coeff05sm=abs(smooth(v1coeff05,5))
w01=where(v1coeff01sm gt v1coeff03sm and v1coeff01sm gt v1coeff05sm)
w03=where(v1coeff03sm gt v1coeff01sm and v1coeff03sm gt v1coeff05sm)
w05=where(v1coeff05sm gt v1coeff01sm and v1coeff05sm gt v1coeff03sm)

c01gtarrv1=fltarr(ntime)*1./0. & c01gtarrv1[w01]=plevel
c03gtarrv1=fltarr(ntime)*1./0. & c03gtarrv1[w03]=plevel
c05gtarrv1=fltarr(ntime)*1./0. & c05gtarrv1[w05]=plevel









yr=(abtimarr-anytim('01-jan-1996 00:00:00'))/3600./24./365.


;VERSION 1 - PLOT 6. - plot polar fields and dipole sphere harmonics---------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_6_polarb_sphereharm.eps'
setcolors,/sys,/silent,/quiet
;mintim=min(tflux27dseries)
mintim=min(ABTIMARR)

linecol=0

;polar plot
setcolors,/sys
utplot,timrebin-mintim,npole6065,mintim,yran=[-9.99,10], pos=[.12,.66,.95,.95],ps=10,xran=minmax(timrebin)-mintim,/xsty,ytit='Polar Fields [G]',xtit='',xtickname=strarr(10)+' '
oplot,timrebin-mintim,npole6065,color=!red,ps=10
oplot,timrebin-mintim,npole6570,color=!red,ps=10,lines=2
;oplot,timrebin-mintim,npole7075,color=!red,ps=10
oplot,timrebin-mintim,spole6065,color=!blue,ps=10
oplot,timrebin-mintim,spole6570,color=!blue,ps=10,lines=2
;oplot,timrebin-mintim,spole7075,color=!blue,ps=10
hline,0,color=!gray

legend,['Lat. 60-65','Lat. 65-70'],/top,/right,color=[linecol,linecol],lines=[0,2]


;spherical harmonic plot
;monopole
utplot,ABTIMARR-mintim,reform(v1COEFF02),mintim,pos=[.12,.37,.95,.66],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-1.99,2],ytit='  V1 Harmonic Coeff. (R=R!dSun!n)',xtit='',xtickname=strarr(10)+' '
;3 bands total
oplot,ABTIMARR-mintim,reform(v1COEFF00),ps=10,color=!gray
;5 bands total
oplot,ABTIMARR-mintim,reform(v1COEFF04),ps=10,lines=2
hline,0,color=!gray

oplot,ABTIMARR-mintim,c02gtarrv2,thick=10
oplot,ABTIMARR-mintim,c04gtarrv2,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevel-0.1,'V2',/data

oplot,ABTIMARR-mintim,c02gtarrv1-0.2,thick=10
oplot,ABTIMARR-mintim,c04gtarrv1-0.2,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevel-0.32,'V1',/data

legend,[textoidl('C_{0,0}'),textoidl('C_{2,0}'),textoidl('C_{4,0}')],/bottom,/right,color=[!gray,linecol,linecol],lines=[0,0,2]

;pure dipole
utplot,ABTIMARR-mintim,reform(v1COEFF01),mintim,pos=[.12,.1,.95,.37],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-2,2],ytit='V1 Harmonic Coeff. (R=R!dSun!n)'
;two bands in each hemisphere
oplot,ABTIMARR-mintim,reform(v1COEFF03),ps=10,color=!gray
;three bands in each hemisphere (AR, decay, pole)
oplot,ABTIMARR-mintim,reform(v1COEFF05),ps=10,lines=2

;oplot,ABTIMARR-mintim,reform(v1COEFF02),ps=10,color=!red ;non hemisherically symmetric
;oplot,ABTIMARR-mintim,reform(v1COEFF04),ps=10,color=!green ;non hemisherically symmetric
hline,0,color=!gray

oplot,ABTIMARR-mintim,c01gtarrv2,thick=10
oplot,ABTIMARR-mintim,c03gtarrv2,thick=10,color=!gray
oplot,ABTIMARR-mintim,c05gtarrv2,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevel-0.1,'V2',/data

oplot,ABTIMARR-mintim,c01gtarrv1-0.2,thick=10
oplot,ABTIMARR-mintim,c03gtarrv1-0.2,thick=10,color=!gray
oplot,ABTIMARR-mintim,c05gtarrv1-0.2,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevel-0.32,'V1',/data

legend,[textoidl('C_{1,0}'),textoidl('C_{3,0}'),textoidl('C_{5,0}')],/bottom,/right,color=[linecol,!gray,linecol],lines=[0,0,2]

closeplotenv

;----------------------------------------------------------------------------->





;yaxis value to over plot and indicate the COEFF of greatest magnitude
plevelr175=0.45

;Smooth by ~3 rotations and determine which is strongest at each time
v1coeff02r175sm=abs(smooth(v1coeff02r175,5))
v1coeff04r175sm=abs(smooth(v1coeff04r175,5))
w02=where(v1coeff02r175sm gt v1coeff04r175sm)
w04=where(v1coeff04r175sm gt v1coeff02r175sm)

c02gtarrv1r175=fltarr(ntime)*1./0. & c02gtarrv1r175[w02]=plevelr175
c04gtarrv1r175=fltarr(ntime)*1./0. & c04gtarrv1r175[w04]=plevelr175

v1coeff01r175sm=abs(smooth(v1coeff01r175,5))
v1coeff03r175sm=abs(smooth(v1coeff03r175,5))
v1coeff05r175sm=abs(smooth(v1coeff05r175,5))
w01=where(v1coeff01r175sm gt v1coeff03r175sm and v1coeff01r175sm gt v1coeff05r175sm)
w03=where(v1coeff03r175sm gt v1coeff01r175sm and v1coeff03r175sm gt v1coeff05r175sm)
w05=where(v1coeff05r175sm gt v1coeff01r175sm and v1coeff05r175sm gt v1coeff03r175sm)

c01gtarrv1r175=fltarr(ntime)*1./0. & c01gtarrv1r175[w01]=plevelr175
c03gtarrv1r175=fltarr(ntime)*1./0. & c03gtarrv1r175[w03]=plevelr175
c05gtarrv1r175=fltarr(ntime)*1./0. & c05gtarrv1r175[w05]=plevelr175





;VERSION 1 - PLOT 6.1 - plot sphere harmonics for R=1.75---------------------------------->

setplotenv,xs=15,ys=15,/ps,file=plotpath+'plot_6_1_polarb_sphereharm_r1_75.eps'
setcolors,/sys,/silent,/quiet
;mintim=min(tflux27dseries)
mintim=min(ABTIMARR)

linecol=0

;spherical harmonic plot
;monopole
utplot,ABTIMARR-mintim,reform(v1coeff02r175),mintim,pos=[.12,.55,.95,.95],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-.5,.5],/ysty,ytit='V1 Harmonic Coeff. (R=1.75R!dSun!n)',xtit='',xtickname=strarr(10)+' '
;3 bands total
oplot,ABTIMARR-mintim,reform(v1coeff00r175),ps=10,color=!gray
;5 bands total
oplot,ABTIMARR-mintim,reform(v1coeff04r175),ps=10,lines=2
hline,0,color=!gray
oplot,ABTIMARR-mintim,c02gtarrv2r175,thick=10
oplot,ABTIMARR-mintim,c04gtarrv2r175,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevelr175-0.025,'V2',/data

oplot,ABTIMARR-mintim,c02gtarrv1r175-0.05,thick=10
oplot,ABTIMARR-mintim,c04gtarrv1r175-0.05,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevelr175-0.075,'V1',/data

legend,[textoidl('C_{0,0}'),textoidl('C_{2,0}'),textoidl('C_{4,0}')],/bottom,/right,color=[!gray,linecol,linecol],lines=[0,0,2]

;pure dipole
utplot,ABTIMARR-mintim,reform(v1coeff01r175),mintim,pos=[.12,.1,.95,.55],ps=10,/noerase,xran=minmax(timrebin)-mintim,/xsty,yran=[-.5,.5],/ysty,ytit='V1 Harmonic Coeff. (R=1.75R!dSun!n)'
;two bands in each hemisphere
oplot,ABTIMARR-mintim,reform(v1coeff03r175),ps=10,color=!gray
;three bands in each hemisphere (AR, decay, pole)
oplot,ABTIMARR-mintim,reform(v1coeff05r175),ps=10,lines=2

;oplot,ABTIMARR-mintim,reform(COEFF02),ps=10,color=!red ;non hemisherically symmetric
;oplot,ABTIMARR-mintim,reform(COEFF04),ps=10,color=!green ;non hemisherically symmetric
hline,0,color=!gray
oplot,ABTIMARR-mintim,c01gtarrv2r175,thick=10
oplot,ABTIMARR-mintim,c03gtarrv2r175,thick=10,color=!gray
oplot,ABTIMARR-mintim,c05gtarrv2r175,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevelr175-0.025,'V2',/data

oplot,ABTIMARR-mintim,c01gtarrv1r175-0.05,thick=10
oplot,ABTIMARR-mintim,c03gtarrv1r175-0.05,thick=10,color=!gray
oplot,ABTIMARR-mintim,c05gtarrv1r175-0.05,thick=5,lines=2
xyouts,anytim('1-jan-2010')-mintim,plevelr175-0.075,'V1',/data

legend,[textoidl('C_{1,0}'),textoidl('C_{3,0}'),textoidl('C_{5,0}')],/bottom,/right,color=[linecol,!gray,linecol],lines=[0,0,2]

closeplotenv

















