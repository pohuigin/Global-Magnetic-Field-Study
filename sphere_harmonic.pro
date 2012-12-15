;From DeRosa PFSS Manual
;THETA=LATITUDE, PHI=LONGITUDE
;------------------------------------------------------------------------->

function sphere_harmonic_coeff, ina_lm, inb_lm, mm=inmm, ll=inll, rr=inrr, use_br=use_br, use_phi=use_phi

a_lm=ina_lm
b_lm=inb_lm
if n_elements(inrr) gt 0 then rr=inrr else rr=1.
mm=inmm
ll=inll

if not keyword_set(use_br) and not keyword_set(use_phi) then use_br=1

;phi_lm=(A_lm[i,j]*rr^(ll[i]) + B_lm[i,j]*rr^(-1.-ll[i]))*Y_lm ;eq. 3
if keyword_set(use_phi) then $
	coeff=A_lm[ll,mm,*]*rr^(ll) + B_lm[ll,mm,*]*rr^(-1.-ll)

;bb_r_lm=(-1.)*Y_lm*(A_lm[i,j]*ll[i]*rr^(ll[i]-1.) - B_lm[i,j]*(ll[i]+1.*rr^(-ll[i]-2.))) ;eq. 13
if keyword_set(use_br) then $
	coeff=(-1.)*(A_lm[ll,mm,*]*ll*rr^(ll-1.) - B_lm[ll,mm,*]*(ll+1.*rr^(-ll-2.)))



return,coeff

end

;------------------------------------------------------------------------->
pro sphere_harmonic_plot

plotp='~/science/papers/active_regions_2_cycle/images/pfss/'

harmfile='~/science/papers/active_regions_2_cycle/data/pfss_monthly_phiab_coeff.sav'
restore,harmfile,/ver
A_lm=PHIATARR
B_lm=PHIBTARR
abtimarr=anytim(date)

;AXISYMMETRIC MODES

;coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1., /use_phi)
;coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1., /use_phi)
;coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1., /use_phi)
;coeff03=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=1., /use_phi)
;coeff04=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=1., /use_phi)
;coeff05=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=1., /use_phi)
;
;save,abtimarr,coeff00,coeff01,coeff02,coeff03,coeff04,coeff05,file='~/science/papers/active_regions_2_cycle/data/pfss_monthly_phiab_plot.sav'

;coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=2.5, /use_phi)
;coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=2.5, /use_phi)
;coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=2.5, /use_phi)
;coeff03=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=2.5, /use_phi)
;coeff04=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=2.5, /use_phi)
;coeff05=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=2.5, /use_phi)

;save,abtimarr,coeff00,coeff01,coeff02,coeff03,coeff04,coeff05,file='~/science/papers/active_regions_2_cycle/data/pfss_monthly_phiab_plot_r2_5.sav'

coeff00=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=0., rr=1.75, /use_phi)
coeff01=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=1., rr=1.75, /use_phi)
coeff02=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=2., rr=1.75, /use_phi)
coeff03=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=3., rr=1.75, /use_phi)
coeff04=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=4., rr=1.75, /use_phi)
coeff05=sphere_harmonic_coeff(a_lm, b_lm, mm=0., ll=5., rr=1.75, /use_phi)

save,abtimarr,coeff00,coeff01,coeff02,coeff03,coeff04,coeff05,file='~/science/papers/active_regions_2_cycle/data/pfss_monthly_phiab_plot_r1_75.sav'


stop

setcolors,/sys,/silent,/quiet
utplot,anytim(nowarr)-min(anytim(nowarr)),coeff00,min(anytim(nowarr)),ps=-1,yrange=[-2,1.5]
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff01,ps=-2,color=!red
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff02,ps=-3,color=!orange
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff03,ps=-4,color=!yellow
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff04,ps=-5,color=!green
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff05,ps=-6,color=!cyan

window_capture,file=plotp+'sphere_harm_axisymmetric_meq0',/png

;NON-AXISYMMETRIC MODES (split longitudinally, m=1)

coeff10=sphere_harmonic_coeff(a_lm, b_lm, mm=1., ll=0., rr=1., /use_phi)
coeff11=sphere_harmonic_coeff(a_lm, b_lm, mm=1., ll=1., rr=1., /use_phi)
coeff12=sphere_harmonic_coeff(a_lm, b_lm, mm=1., ll=2., rr=1., /use_phi)
coeff13=sphere_harmonic_coeff(a_lm, b_lm, mm=1., ll=3., rr=1., /use_phi)
coeff14=sphere_harmonic_coeff(a_lm, b_lm, mm=1., ll=4., rr=1., /use_phi)
coeff15=sphere_harmonic_coeff(a_lm, b_lm, mm=1., ll=5., rr=1., /use_phi)

setcolors,/sys,/silent,/quiet
utplot,anytim(nowarr)-min(anytim(nowarr)),coeff10,min(anytim(nowarr)),ps=-1,yrange=[-2,1.5]
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff11,ps=-2,color=!red
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff12,ps=-3,color=!orange
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff13,ps=-4,color=!yellow
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff14,ps=-5,color=!green
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff15,ps=-6,color=!cyan

window_capture,file=plotp+'sphere_harm_nonaxisymmetric_meq1',/png

;NON-AXISYMMETRIC MODES (split longitudinally, m=2)

coeff20=sphere_harmonic_coeff(a_lm, b_lm, mm=2., ll=0., rr=1., /use_phi)
coeff21=sphere_harmonic_coeff(a_lm, b_lm, mm=2., ll=1., rr=1., /use_phi)
coeff22=sphere_harmonic_coeff(a_lm, b_lm, mm=2., ll=2., rr=1., /use_phi)
coeff23=sphere_harmonic_coeff(a_lm, b_lm, mm=2., ll=3., rr=1., /use_phi)
coeff24=sphere_harmonic_coeff(a_lm, b_lm, mm=2., ll=4., rr=1., /use_phi)
coeff25=sphere_harmonic_coeff(a_lm, b_lm, mm=2., ll=5., rr=1., /use_phi)

setcolors,/sys,/silent,/quiet
utplot,anytim(nowarr)-min(anytim(nowarr)),coeff20,min(anytim(nowarr)),ps=-1,yrange=[-2,1.5]
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff21,ps=-2,color=!red
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff22,ps=-3,color=!orange
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff23,ps=-4,color=!yellow
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff24,ps=-5,color=!green
oplot,anytim(nowarr)-min(anytim(nowarr)),coeff25,ps=-6,color=!cyan

window_capture,file=plotp+'sphere_harm_nonaxisymmetric_meq2',/png

save,nowarr,coeff00,coeff01,coeff02,coeff03,coeff04,coeff05, $
	coeff10,coeff11,coeff12,coeff13,coeff14,coeff15, $
	coeff20,coeff21,coeff22,coeff23,coeff24,coeff25, $
	file='~/science/papers/active_regions_2_cycle/data/sphere_harm_mothly_timeseries.sav'

stop

end

;------------------------------------------------------------------------->

pro sphere_harmonic

plotp='~/science/papers/active_regions_2_cycle/images/pfss/'
harmfile='~/science/papers/active_regions_2_cycle/data/pfss_yearly_phiab_coeff.sav'

restore,harmfile,/ver

ll=[0.,1.,2.,3.,4.]
;ll=[5.];,1.,2.,3.,4.]
;mm=[0.,1.,2.,3.,4.]
rr=1.
lat=findgen(1,180)/179.*!pi;-!pi/2.
lat=rebin(lat,360,180)
lon=findgen(360,1)/359.*2.*!pi;-!pi/2.
lon=rebin(lon,360,180)

;Temporary until i figure out how to use the PFSS calculated coefficients
A_lm=PHIATARR[*,*,0];1.
B_lm=PHIBTARR[*,*,0];1.

for i=0,n_elements(ll)-1 do begin ;loop over l
	;mm can range -ll to ll, but here, just take positive...
	if i eq 0 then mm=0. else mm=findgen(i+1)
	
	j=0 ;just do axisymmetric modes
;	for j=0,n_elements(mm)-1 do begin ;loop over m
				
		Y_lm=spher_harm(lat,lon,ll[i],mm[j]);C_lm*P_lm*exp(imaginary(mm*lon)) ;eq. 4
		;C_lm=(-1.)^(mm)*sqrt((2.*ll+1.)*factorial(ll-mm)/(4.*!pi)/factorial(ll+mm))
		;P_lm=*cos(lat)
		
		phi_lm=(A_lm[i,j]*rr^(ll[i]) + B_lm[i,j]*rr^(-1.-ll[i]))*Y_lm ;eq. 3
		
		bb=-deriv(phi_lm) ;eq. 1		
;		plot_image,bb,tit='L='+strtrim(i,2)+' M='+strtrim(j,2)
;		stop
		
		bb_r_lm=(-1.)*Y_lm*(A_lm[i,j]*ll[i]*rr^(ll[i]-1.) - B_lm[i,j]*(ll[i]+1.*rr^(-ll[i]-2.))) ;eq. 13

		plot_image,bb_r_lm,tit='L='+strtrim(i,2)+' M='+strtrim(j,2)
		stop
		
		;bb_lat_lm=(-1./(rr*sin(lat))*Y_lm*(R_lm*(ll-1.)*(A_lm[])
		
		;bb_lon_lm=

;plot onto sphere
		MESH_OBJ, 4, vertices, polygons, REPLICATE(0.25, 101, 101), /CLOSED
  	  	oModel = OBJ_NEW('IDLgrModel')  
		oPalette = OBJ_NEW('IDLgrPalette')  
		oPalette -> LOADCT, 0  
		oPalette -> SetRGB, 255, 255, 255, 255  
		oImage = OBJ_NEW('IDLgrImage', bytscl(congrid(bb_r_lm,180,180)), PALETTE = oPalette)
		vector = FINDGEN(101)/100.  
		texure_coordinates = FLTARR(2, 101, 101)  
		texure_coordinates[0, *, *] = vector # REPLICATE(1., 101)  
		texure_coordinates[1, *, *] = REPLICATE(1., 101) # vector
		oPolygons = OBJ_NEW('IDLgrPolygon', SHADING = 0, $  
			DATA = vertices, POLYGONS = polygons, $  
			COLOR = [255, 255, 255], $  
			TEXTURE_COORD = texure_coordinates, $  
			TEXTURE_MAP = oImage, /TEXTURE_INTERP)
		oModel -> ADD, oPolygons  
		oModel -> ROTATE, [1, 0, 0], -90  
		oModel -> ROTATE, [0, 1, 0], -90
		XOBJVIEW, oModel, /BLOCK 	
		
		window_Capture,file=plotp+'sphere_harm_vis_l'+strtrim(i,2)+'m'+strtrim(j,2),/png

stop
		
		OBJ_DESTROY, [oModel, oImage, oPalette] 
	
;	endfor
endfor






end