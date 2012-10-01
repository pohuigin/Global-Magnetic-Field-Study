pro smart_cycle_modulation

datapath='~/science/data/cycle23_sav/smart_paper_2/'
restore,datapath+'butterfly_diags_lat.sav'
!x.margin=[10,2]
!y.margin=[3.5,.5]

xbin=14./365.

;MODULATION R-VALUE IMAGE North----------------------------------------->
fxarr=findgen(n_elements(rvalimage[*,0]))*xbin*365.
fyarr=smooth(total(rvalimage[*,90:179],2),4)

;frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years',ytit='Flux [Mx]'


period=1./xpeak[[9,6]]
xyouts,.6,.95,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.92,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.89,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.86,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.83,strjoin([strtrim(14.*2./(30.4167),2), 'Nyq in months'],' '),/norm
;xyouts,.6,.83,strjoin([strtrim(14.*2.,2),strtrim(14.*2./365.,2),strtrim(14.*2./(30.4167),2),'nyquist in days years months'],2),/norm
window_capture,file=datapath+'dft_rval_north',/png
stop

;MODULATION R-VALUE IMAGE South----------------------------------------->
fxarr=findgen(n_elements(rvalimage[*,0]))*xbin*365.
fyarr=smooth(total(rvalimage[*,0:89],2),4)

frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years'

period=1./xpeak[[11,10,9]]
xyouts,.6,.9,strjoin(['Frequesncy1','Frequency2','Frequency3'],' '),/norm
xyouts,.6,.87,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.84,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.81,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.78,strjoin([strtrim(14.*2./(30.4167),2),'nyquist in months'],' '),/norm
window_capture,file=datapath+'dft_rval_south',/png
stop

;MODULATION WLSG IMAGE North----------------------------------------->
fxarr=findgen(n_elements(wlsgimage[*,0]))*xbin*365.
fyarr=smooth(total(wlsgimage[*,90:179],2),4)

;frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years',ytit='Flux [Mx]'

period=1./xpeak[[7,10]]
xyouts,.6,.95,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.92,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.89,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.86,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.83,strjoin([strtrim(14.*2./(30.4167),2), 'Nyq in months'],' '),/norm
;xyouts,.6,.83,strjoin([strtrim(14.*2.,2),strtrim(14.*2./365.,2),strtrim(14.*2./(30.4167),2),'nyquist in days years months'],2),/norm
window_capture,file=datapath+'dft_wlsg_north',/png
stop

;MODULATION WLSG IMAGE South----------------------------------------->
fxarr=findgen(n_elements(wlsgimage[*,0]))*xbin*365.
fyarr=smooth(total(wlsgimage[*,0:89],2),4)

frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years'

period=1./xpeak[[12,13,14,16]]
xyouts,.6,.9,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.87,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.84,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.81,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.78,strjoin([strtrim(14.*2./(30.4167),2),'nyquist in months'],' '),/norm
window_capture,file=datapath+'dft_wlsg_south',/png
stop

;MODULATION FLUX IMAGE North----------------------------------------->
fxarr=findgen(n_elements(fluximage[*,0]))*xbin*365.
fyarr=smooth(total(fluximage[*,90:179],2),4)

;frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years',ytit='Flux [Mx]'

period=1./xpeak[[6,7,11,16]]
xyouts,.6,.95,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.92,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.89,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.86,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.83,strjoin([strtrim(14.*2./(30.4167),2), 'Nyq in months'],' '),/norm
;xyouts,.6,.83,strjoin([strtrim(14.*2.,2),strtrim(14.*2./365.,2),strtrim(14.*2./(30.4167),2),'nyquist in days years months'],2),/norm
window_capture,file=datapath+'dft_flux_north',/png
stop

;MODULATION FLUX IMAGE South----------------------------------------->
fxarr=findgen(n_elements(fluximage[*,0]))*xbin*365.
fyarr=smooth(total(fluximage[*,0:89],2),4)

frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years'

period=1./xpeak[[13,15]]
xyouts,.6,.9,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.87,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.84,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.81,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.78,strjoin([strtrim(14.*2./(30.4167),2),'nyquist in months'],' '),/norm
window_capture,file=datapath+'dft_flux_south',/png
stop

;MODULATION FLARE INDEX IMAGE North----------------------------------------->
fxarr=findgen(n_elements(flrindeximage[*,0]))*xbin*365.
fyarr=smooth(total(flrindeximage[*,90:179],2),4)

;frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years',ytit='Flux [Mx]'

period=1./xpeak[[3,8,10,12]]
xyouts,.6,.95,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.92,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.89,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.86,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.83,strjoin([strtrim(14.*2./(30.4167),2), 'Nyq in months'],' '),/norm
;xyouts,.6,.83,strjoin([strtrim(14.*2.,2),strtrim(14.*2./365.,2),strtrim(14.*2./(30.4167),2),'nyquist in days years months'],2),/norm
window_capture,file=datapath+'dft_flrindex_north',/png
stop

;MODULATION FLARE INDEX IMAGE South----------------------------------------->
fxarr=findgen(n_elements(flrindeximage[*,0]))*xbin*365.
fyarr=smooth(total(flrindeximage[*,0:89],2),4)

frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years'

period=1./xpeak[[11,13]]
xyouts,.6,.9,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.87,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.84,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.81,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.78,strjoin([strtrim(14.*2./(30.4167),2),'nyquist in months'],' '),/norm
window_capture,file=datapath+'dft_flrindex_south',/png
stop

;MODULATION FLARE MX IMAGE North----------------------------------------->
fxarr=findgen(n_elements(flrmximage[*,0]))*xbin*365.
fyarr=smooth(total(flrmximage[*,90:179],2),4)

;frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years',ytit='Flux [Mx]'

period=1./xpeak[[15,16]]
xyouts,.6,.95,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.92,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.89,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.86,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.83,strjoin([strtrim(14.*2./(30.4167),2), 'Nyq in months'],' '),/norm
;xyouts,.6,.83,strjoin([strtrim(14.*2.,2),strtrim(14.*2./365.,2),strtrim(14.*2./(30.4167),2),'nyquist in days years months'],2),/norm
window_capture,file=datapath+'dft_flrmx_north',/png
stop

;MODULATION FLARE MX IMAGE South----------------------------------------->
fxarr=findgen(n_elements(flrmximage[*,0]))*xbin*365.
fyarr=smooth(total(flrmximage[*,0:89],2),4)

frange=[0,0.00555556*2.] ;0 dayHz, 1 revolution/year in dayHz
dft_wrap, fxarr, fyarr, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years'

period=1./xpeak[[12,16]]
xyouts,.6,.9,strjoin(['Frequesncy1','Frequency2'],' '),/norm
xyouts,.6,.87,strjoin([strtrim(period,2), 'days'],' '),/norm
xyouts,.6,.84,strjoin([strtrim(period/365.,2), 'years'],' '),/norm
xyouts,.6,.81,strjoin([strtrim(period/(30.4167),2),'months'],' '),/norm
xyouts,.6,.78,strjoin([strtrim(14.*2./(30.4167),2),'nyquist in months'],' '),/norm
window_capture,file=datapath+'dft_flrmx_south',/png
stop


restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday.sav',/ver
restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_tim.sav',/ver
restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_carrarr.sav',/ver

;1 DAY BINNED FLUX MODULATION
smart_bin_time, arstr_arr.bflux, tim1pd, odata, otime, bin=1./365.,/total
otime=findgen(n_elements(odata))*(1./365.)
plot,(real_part(fft(odatasmth))^2.)[0:2373],xran=[10,2500]
dft_wrap, otime, odata, xout, yout, frange=frange, /plot, wpeak=wpeak, xpeak=xpeak,xtit='Years',ytit='Flux [Mx]'

stop



;Correlation between flares and flux
;FLARE INDEX
print,correlate(smooth(total(flrindeximage,2),4),smooth(total(rvalimage,2),4));     0.815566
print,correlate(smooth(total(flrindeximage,2),4),smooth(total(wlsgimage,2),4));     0.768988
print,correlate(smooth(total(flrindeximage,2),4),smooth(total(fluximage,2),4));     0.714068
print,correlate(smooth(total(flrindeximage,2),4),smooth(total(areaimage,2),4));     0.667949
;FLARE IMAGE
print,correlate(smooth(total(flrmximage,2),4),smooth(total(wlsgimage,2),4));     0.880441
print,correlate(smooth(total(flrmximage,2),4),smooth(total(fluximage,2),4));     0.874823
print,correlate(smooth(total(flrmximage,2),4),smooth(total(areaimage,2),4));     0.872659






end