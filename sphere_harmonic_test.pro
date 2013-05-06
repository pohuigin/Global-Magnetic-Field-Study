;routine to test plot the v1 and v2 PFSS runs and compare the a,b, and C's...
;@sphere_harmonic_test.pro

.r sphere_harmonic

;VERSION 2-------------->

harmfile='~/science/papers/active_regions_2_cycle/data/pfss_monthly_phiab_coeff_v2.sav'
restore,harmfile,/ver
A_lm=PHIATARR
B_lm=PHIBTARR
abtimarr=anytim(date)

v2coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1, /use_phi,outa=a00, outb=b00)
v2coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1, /use_phi,outa=a01, outb=b01)
v2coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1, /use_phi,outa=a02, outb=b02)

yr=(abtimarr-anytim('01-jan-1996 00:00:00'))/3600./24./365.

window,xs=800,ys=800
!p.multi=[0,1,3]
!y.margin=[3,0]
!p.charsize=2
!y.range=[-2,2]
!x.style=1

plot,yr,a00,ytit='m0,l0',ps=10
oplot,yr,b00,lines=2,ps=10                
plot,yr,a01,ytit='m0,l1',ps=10
oplot,yr,b01,lines=2,ps=10                
plot,yr,a02,ytit='m0,l2',ps=10
oplot,yr,b02,lines=2,ps=10                

imgpath='~/science/papers/active_regions_2_cycle/images/pfss_test_v1v2/'
window_capture,file=imgpath+'test_v2_ab'

stop

;VERSION 1------------------>

harmfile='~/science/papers/active_regions_2_cycle/data/pfss_monthly_phiab_coeff.sav'
restore,harmfile,/ver
A_lm=PHIATARR
B_lm=PHIBTARR
abtimarr=anytim(date)

v1coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1, /use_phi,outa=a00, outb=b00)
v1coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1, /use_phi,outa=a01, outb=b01)
v1coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1, /use_phi,outa=a02, outb=b02)

yr=(abtimarr-anytim('01-jan-1996 00:00:00'))/3600./24./365.

window,xs=800,ys=800
!p.multi=[0,1,3]
!y.margin=[3,0]
!p.charsize=2
!y.range=[-2,2]
!x.style=1

plot,yr,a00,ytit='m0,l0',ps=10
oplot,yr,b00,lines=2,ps=10                
plot,yr,a01,ytit='m0,l1',ps=10
oplot,yr,b01,lines=2,ps=10                
plot,yr,a02,ytit='m0,l2',ps=10
oplot,yr,b02,lines=2,ps=10                

imgpath='~/science/papers/active_regions_2_cycle/images/pfss_test_v1v2/'
window_capture,file=imgpath+'test_v1_ab'

stop

;COMPARE V1 and V2

!p.multi=[0,1,2]
!y.range=[-2,2]
!x.style=1
plot,yr,v1coeff01,ytit='VERSION 1',ps=10
oplot,yr,v1coeff00,ps=10,color=180
oplot,yr,v1coeff02,lines=2,ps=10                
hline,0,color=180
plot,yr,v2coeff01,ytit='VERSION 2',ps=10
oplot,yr,v2coeff00,ps=10,color=180
oplot,yr,v2coeff02,lines=2,ps=10  
hline,0,color=180

imgpath='~/science/papers/active_regions_2_cycle/images/pfss_test_v1v2/'
window_capture,file=imgpath+'test_v1_vs_v2'


