pro run_smart_bin2d_2, res0tore=res0tore, res1tore=res1tore

;restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday.sav',/ver
;restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_tim.sav',/ver ;checked! in units of years
;restore,'~/science/data/smart_sf_smart_compare/arstr_arr_1pday_carrarr.sav',/ver

restore,'~/science/data/goes_event_list/goes_flare_struct.sav',/ver

if not keyword_set(res0tore) then begin
restore,'~/science/papers/active_regions_2_cycle/images/smart_sav_all_20101018_0836.sav',/ver
dates=datearr('19960101','20080101')
ardates=time2file(anytim(arstruct_all.time,/vms),/date)
arstruct_arr=arstruct_all[0]
ardates_arr=''
for i=0l,n_elements(dates)-1 do begin
	wdate=where(ardates eq dates[i])
	if wdate[0] ne -1 then begin
		thisartime=arstruct_all[wdate].time 
		wfirst=where(thisartime eq thisartime[uniq(thisartime)])
		arstruct_arr=[arstruct_arr,(arstruct_all[wdate])]
		ardates_arr=[ardates_arr,(arstruct_all[wdate[0]].time)[wfirst]]
	endif
endfor
arstruct_arr=arstruct_arr[1:*]
ardates_arr=ardates_arr[1:*]
artims_arr=anytim(ardates_arr)
save,artims_arr,ardates_arr,arstruct_arr,file='~/science/papers/active_regions_2_cycle/images/smart_sav_1pday.sav',/ver
endif else restore,'~/science/papers/active_regions_2_cycle/images/smart_sav_1pday.sav',/ver

;to be subtrackted from all time arrays:
mintim=0.

;check carringlon longitude calculation
;carlonarr2=tim2carr(arstr_arr.time)+arstr_arr.hglon
;wclt0=where(carlonarr2 lt 0)
;if wclt0[0] ne -1 then carlonarr2[wclt0]=carlonarr2[wclt0]+360.
;wcge360=where(carlonarr2 ge 360)
;if wcge360[0] ne -1 then carlonarr2[wcge360]=carlonarr2[wcge360]-360.
;stop

datapath='~/science/data/cycle23_sav/smart_paper_2/'
paper2path='~/science/papers/active_regions_2_cycle/images/'

goto,skipto2wkbin

;-------------------------------------------------------------------------------------------------->
;LONGITUDE REBIN TO 2weeks------------------------------------------------------------------------->
;xbin=.08333 ;month bin
xbin=.0383562 ;2 week bin = 14./365.
ybin=1.
xran=[0,11.2]
yran=[0,360]

if not keyword_set(res1tore) then begin 
;Run longitude images
xx=smart_bin2d(xbin=xbin,ybin=ybin,xran=xran,yran=yran,/years,struct=arstruct_arr, $
	flrstruct=astrpos,flrtimarr=FLRPOSSTR.tim/3600./24./365.-mintim, flrposstr=FLRPOSSTR, $
	/longitude,longarr=CARLONARR,timarr=tim1pd,path=datapath,outfile=outfile)
restore,outfile,/ver
endif else restore,datapath+'butterfly_diags_lon.sav',/ver

setplotenv,/xwin & !p.color=0 & !p.background=255
window,8,xs=1000,ys=700

;PLOT carlon flares complexity----------------------------------------->
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

;IMAGE carlon flares complexity----------------------------------------->
loadct,0
plot_image,(rvalimage)^.2,/nosq,origin=[0,0],scale=[xbin,ybin],xtit='years since 19970101',ytit='carrington longitude',tit='R value' 
plotsym,0,.5,/fill
setcolors,/sys
window_capture,file=datapath+'carlon_rvalue_mxflares_image0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,FLRPOSSTR[wmx].carlon,ps=8,color=!red
window_capture,file=datapath+'carlon_rvalue_mxflares_image1',/png

loadct,0
plot_image,(wlsgimage)^.2,/nosq,origin=[0,0],scale=[xbin,ybin],xtit='years since 19970101',ytit='carrington longitude',tit='WLSG' 
plotsym,0,.5,/fill
setcolors,/sys
window_capture,file=datapath+'carlon_wlsg_mxflares_image0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,FLRPOSSTR[wmx].carlon,ps=8,color=!red
window_capture,file=datapath+'carlon_wlsg_mxflares_image1',/png
stop

skipto2wkbin:
;-------------------------------------------------------------------------------------------------->
;LATITUDE REBIN TO 2weeks-------------------------------------------------------------------------->
xbin=.0383562 ;2 week bin = 14./365.
ybin=1.
xran=[anytim('1-jan-1997'),anytim('1-jan-2009')]/(3600.*24.*365.)
yran=[-90,90]
if not keyword_set(res1tore) then begin
;Run latitude images
;xx=smart_bin2d(xbin=xbin,ybin=ybin,xran=xran,yran=yran,/years,struct=arstr_arr, $
;	flrstruct=astrpos,flrtimarr=FLRPOSSTR.tim/3600./24./365.-mintim, flrposstr=FLRPOSSTR, $
;	/latitude,timarr=tim1pd,path=datapath,outfile=outfile)
tim1pd=anytim(arstruct_arr.time)
w2wkbin=where(tim1pd gt anytim('1-jan-1997') and tim1pd lt anytim('1-jan-2009')+11.2*365.*24.*3600.)
;tim1pd=(tim1pd-min(tim1pd))/3600.*24.*365.
tim1pd=(tim1pd)/(3600.*24.*365.)
xx=smart_bin2d(xbin=xbin,ybin=ybin,xran=xran,yran=yran,/years,struct=arstruct_arr, $
	flrstruct=astrpos,flrtimarr=FLRPOSSTR.tim/3600./24./365.-mintim, flrposstr=FLRPOSSTR, $
	/latitude,timarr=tim1pd,path=datapath,outfile=outfile)
restore,outfile,/ver
endif else restore,datapath+'butterfly_diags_lat.sav'

;MAKE PLOTs
window,xs=1000,ys=700
!p.multi=[0,1,2]

;flux balance avged (MORE + / MORE -) in all ARs thresholded to binary
setplotenv,file=paper2path+'butterfly_2wk_polarity_only.eps',/ps,xs=15,ys=10
!p.multi=0
loadct,0
;plot_image,(rebin(FLUXBALIMAGE,292.*2.,360,/samp)) < .001> (-.001),/nosq
plot_image,(FLUXBALIMAGE[*,30:149]) < .001> (-.001),/nosq,origin=[0,-60],scale=[xbin,1.]
setcolors,/sys
legend,['More (+)','More (-)'],psym=[4,4],colors=[!white,!black],/top,/right
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_only.eps '+paper2path+'butterfly_2wk_polarity_only.pdf'

setplotenv,file=paper2path+'butterfly_2wk_polarity_marg_only.eps',/ps,xs=10,ys=6
setcolors,/sys
plot,findgen(180)-90.,abs(total(FLUXBALIMAGE<0,1))/292.,xran=[-60,60],/xsty,ytit='<'+textoidl('\Phi')+' Imb.>',xtickname=strarr(10)+' '
oplot,findgen(180)-90.,abs(total(FLUXBALIMAGE>0,1))/292.,color=!red
legend,['More (+)','More (-)'],psym=[4,4],colors=[!red,!black],/top,/right
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_marg_only.eps '+paper2path+'butterfly_2wk_polarity_marg_only.pdf'

;flux balance avged (MORE + / MORE -) in all ARs
setplotenv,file=paper2path+'butterfly_2wk_polarity_ammt.eps',/ps,xs=15,ys=10
!p.multi=0
loadct,0
;plot_image,(rebin(FLUXBALIMAGE,292.*2.,360,/samp)) < .001> (-.001),/nosq
plot_image,(FLUXBALIMAGE[*,30:149]),/nosq,origin=[0,-60],scale=[xbin,1.]
setcolors,/sys
legend,['More (+)','More (-)'],psym=[4,4],colors=[!white,!black],/top,/right
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_ammt.eps '+paper2path+'butterfly_2wk_polarity_ammt.pdf'

;NET flux balance in all ARs with time series (diff of tot phi+ and tot phi-)
setplotenv,file=paper2path+'butterfly_2wk_polarity_net.eps',/ps,xs=15,ys=10
!p.multi=[0,1,2]
loadct,0 & tvlct,rr0,gg0,bb0,/get
plot_image,(fluxposimage-fluxnegimage)[*,30:149] < (1d6) > (-1d6),/nosq,origin=[0,-60],scale=[xbin,1.],ytit='Latitude [deg]',xtickname=strarr(10)+' ',ymarg=[0,1]
hline,0,color=0
setcolors,/sys
plot,xarr+1979.,smooth(total((fluxposimage-fluxnegimage)[*,90:179],2),5)*1d16,ps=10,yran=[-3d24,3d24],thick=5,/xsty,/ysty,ymarg=[5,0],xtit='Years',ytit=textoidl('(\Phi_{+}+\Phi_{-})')+' [Mx]'
hline,0,color=0
oplot,xarr+1979.,smooth(total((fluxposimage-fluxnegimage)[*,0:89],2),5)*1d16,ps=10,color=!red,thick=5,lines=2
legend,['North','South'],lines=[0,2],colors=[!black,!red],thick=[5,5],/top,/right,/fill
color_table, minmax(((fluxposimage-fluxnegimage)[*,30:149])*1d16 < (1d23) > (-1d23)),[.7,.95],[.9,.95],rr=rr0,gg=gg0,bb=bb0,title=textoidl('(\Phi_{+}+\Phi_{-})')+' [Mx]'
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_net.eps '+paper2path+'butterfly_2wk_polarity_net.pdf'

;no butterfly
!p.multi=0
setplotenv,file=paper2path+'butterfly_2wk_polarity_net_nobutt.eps',/ps,xs=15,ys=6
plot,xarr+1979.,smooth(total((fluxposimage-fluxnegimage)[*,90:179],2),5)*1d16,ps=10,yran=[-3d24,3d24],thick=5,/xsty,xran=[1996,2009],/ysty,ymarg=[5,0],xtit='Years',ytit=textoidl('(\Phi_{+}+\Phi_{-})')+' [Mx]',ytickname=strarr(10)+' '
hline,0,color=0
oplot,xarr+1979.,smooth(total((fluxposimage-fluxnegimage)[*,0:89],2),5)*1d16,ps=10,color=!red,thick=5,lines=2
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_net_nobutt.eps '+paper2path+'butterfly_2wk_polarity_net_nobutt.pdf'


stop

;NET flux balance in unipole ARs (diff of tot phi+ and tot phi-) and d(phi)/dt AVGED (N emerge/ N decay) !!!
setplotenv,file=paper2path+'butterfly_2wk_polarity_net_uni_dphidt_avg.eps',/ps,xs=15,ys=10
!p.multi=[0,1,2]
loadct,0 & tvlct,rr0,gg0,bb0,/get
plot_image,(fluxunipolpos-fluxunipolneg)[*,30:149] < (1d5) > (-1d5),/nosq,origin=[0,-60],scale=[xbin,1.],ytit='Latitude [deg]',xtickname=strarr(10)+' ',ymarg=[2,1]
hline,0,color=0
plot_image,(fluxuniemergeavg)[*,30:149] < (1) > (-1),/nosq,origin=[0,-60],scale=[xbin,1.],ytit='Latitude [deg]',xtit='Years since 1997',ymarg=[3,-2]
hline,0,color=0
color_table, minmax(((fluxunipolpos-fluxunipolneg)[*,30:149])*1d16 < (1d23) > (-1d23)),[.7,.95],[.9,.95],rr=rr0,gg=gg0,bb=bb0,title=textoidl('(\Phi_{+}+\Phi_{-})')+' [Mx]'
color_table, minmax(fluxuniemergeavg[*,30:149]),[.7,.95],[.42,.47],rr=rr0,gg=gg0,bb=bb0,title=textoidl('(N_{emerge}-N_{decay})/N_{tot}')
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_net_uni_dphidt_avg.eps '+paper2path+'butterfly_2wk_polarity_net_uni_dphidt_avg.pdf'

stop

;NET flux balance in unipole ARs (diff of tot phi+ and tot phi-) and d(phi)/dt
setplotenv,file=paper2path+'butterfly_2wk_polarity_net_uni_dphidt.eps',/ps,xs=15,ys=10
!p.multi=[0,1,2]
loadct,0 & tvlct,rr0,gg0,bb0,/get
plot_image,(fluxunipolpos-fluxunipolneg)[*,30:149] < (1d5) > (-1d5),/nosq,origin=[0,-60],scale=[xbin,1.],ytit='Latitude [deg]',xtickname=strarr(10)+' ',ymarg=[2,1]
hline,0,color=0
plot_image,(fluxuniemerge)[*,30:149] < (1) > (-1),/nosq,origin=[0,-60],scale=[xbin,1.],ytit='Latitude [deg]',xtit='Years since 1997',ymarg=[3,-2]
hline,0,color=0
color_table, minmax(((fluxunipolpos-fluxunipolneg)[*,30:149])*1d16 < (1d23) > (-1d23)),[.7,.95],[.9,.95],rr=rr0,gg=gg0,bb=bb0,title=textoidl('(\Phi_{+}+\Phi_{-})')+' [Mx]'
color_table, minmax(((fluxuniemerge)[*,30:149])*1d16 < (1d23) > (-1d23)),[.7,.95],[.42,.47],rr=rr0,gg=gg0,bb=bb0,title=textoidl('(\partial\Phi_{+}/\partial t)')+' [Mx/s]'
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_net_uni_dphidt.eps '+paper2path+'butterfly_2wk_polarity_net_uni_dphidt.pdf'

stop

;NET flux balance in all ARs with time series at HI and LOW lats (diff of tot phi+ and tot phi-)
setplotenv,file=paper2path+'butterfly_2wk_polarity_net_hilo.eps',/ps,xs=15,ys=10
!p.multi=[0,1,2]
setcolors,/sys
plot,xarr,smooth(total((fluxposimage-fluxnegimage)[*,110:179],2),5)*1d16,ps=10,yran=[-1.5d23,1.5d23],thick=5,/xsty,/ysty,ymarg=[0,5],ytit=textoidl('High Lat. (\Phi_{+}+\Phi_{-})')+' [Mx]', xtickname=strarr(10)+' '
hline,0,color=0
oplot,xarr,smooth(total((fluxposimage-fluxnegimage)[*,0:69],2),5)*1d16,ps=10,color=!red,thick=5,lines=2
xyouts,.6,.55,'North: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,110:179])*1d16,form='(E9.2)')+' Mx',/norm
xyouts,.6,.52,'South: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,0:69])*1d16,form='(E9.2)')+' Mx',/norm,color=!red
legend,['North','South'],lines=[0,2],colors=[!black,!red],thick=[5,5],/top,/right
plot,xarr,smooth(total((fluxposimage-fluxnegimage)[*,90:109],2),5)*1d16,ps=10,yran=[-3d23,3d23],thick=5,/xsty,/ysty,ymarg=[5,0],xtit='Years Since 1997',ytit=textoidl('Low Lat. (\Phi_{+}+\Phi_{-})')+' [Mx]'
hline,0,color=0
oplot,xarr,smooth(total((fluxposimage-fluxnegimage)[*,70:89],2),5)*1d16,ps=10,color=!red,thick=5,lines=2
xyouts,.6,.18,'North: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,90:109])*1d16,form='(E9.2)')+' Mx',/norm
xyouts,.6,.15,'South: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,70:89])*1d16,form='(E9.2)')+' Mx',/norm,color=!red
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_net_hilo.eps '+paper2path+'butterfly_2wk_polarity_net_hilo.pdf'

;Cummulative net flux imbalance for high latitudes
nfluxnet=total((fluxposimage-fluxnegimage)[*,110:179],2)*1d16
sfluxnet=total((fluxposimage-fluxnegimage)[*,0:69],2)*1d16
nfluxcumm=fltarr(n_elements(xarr)) & nfluxcumm[0]=nfluxnet[0]
sfluxcumm=nfluxcumm & sfluxcumm[0]=nfluxnet[0]
for i=1,n_elements(xarr)-1 do begin
	nfluxcumm[i]=nfluxcumm[i-1]+nfluxnet[i]
	sfluxcumm[i]=sfluxcumm[i-1]+sfluxnet[i]
endfor
setplotenv,file=paper2path+'butterfly_2wk_polarity_cumm_hilo.eps',/ps,xs=15,ys=10
!p.multi=0
setcolors,/sys
plot,xarr,nfluxcumm,ps=10,thick=5,/xsty,ytit=textoidl('Cummulative (\Phi_{+}+\Phi_{-})')+' [Mx]'
hline,0,color=0
oplot,xarr,sfluxcumm,ps=10,color=!red,thick=5,lines=2
;xyouts,.6,.55,'North: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,110:179])*1d16,form='(E9.2)')+' Mx',/norm
;xyouts,.6,.52,'South: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,0:69])*1d16,form='(E9.2)')+' Mx',/norm,color=!red
legend,['North','South'],lines=[0,2],colors=[!black,!red],thick=[5,5],/top,/right
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_cumm_hilo.eps '+paper2path+'butterfly_2wk_polarity_cumm_hilo.pdf'

;Cummulative net UNIPOLE flux imbalance for high latitudes
nfluxnet=total((fluxunipolpos-fluxunipolneg)[*,110:179],2)*1d16
sfluxnet=total((fluxunipolpos-fluxunipolneg)[*,0:69],2)*1d16
nfluxcumm=fltarr(n_elements(xarr)) & nfluxcumm[0]=nfluxnet[0]
sfluxcumm=nfluxcumm & sfluxcumm[0]=nfluxnet[0]
for i=1,n_elements(xarr)-1 do begin
	nfluxcumm[i]=nfluxcumm[i-1]+nfluxnet[i]
	sfluxcumm[i]=sfluxcumm[i-1]+sfluxnet[i]
endfor
setplotenv,file=paper2path+'butterfly_2wk_polarity_cumm_uni_hilo.eps',/ps,xs=15,ys=10
!p.multi=[0,1,2]
setcolors,/sys
plot,xarr,nfluxnet,ps=10,thick=5,/xsty,yran=[-10d22,20d22],ytit=textoidl('Unipolar (\Phi_{+}+\Phi_{-})')+' [Mx]',ymarg=[2,1]
hline,0,color=0
oplot,xarr,sfluxnet,ps=10,color=!red,thick=5,lines=2
plot,xarr,nfluxcumm,ps=10,thick=5,/xsty,yran=[-10d22,20d22],ytit=textoidl('Cummulative Unipolar (\Phi_{+}+\Phi_{-})')+' [Mx]',ymarg=[3,-2]
hline,0,color=0
oplot,xarr,sfluxcumm,ps=10,color=!red,thick=5,lines=2
;xyouts,.6,.55,'North: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,110:179])*1d16,form='(E9.2)')+' Mx',/norm
;xyouts,.6,.52,'South: '+textoidl('(\Phi_{+}+\Phi_{-}) = ')+string(total((fluxposimage-fluxnegimage)[*,0:69])*1d16,form='(E9.2)')+' Mx',/norm,color=!red
legend,['North','South'],lines=[0,2],colors=[!black,!red],thick=[5,5],/top,/right
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_polarity_cumm_uni_hilo.eps '+paper2path+'butterfly_2wk_polarity_cumm_uni_hilo.pdf'


stop

!p.multi=0

;unsigned flux balance avged in all regions
setplotenv,file=paper2path+'butterfly_2wk_fluximb.eps',/ps
loadct,1
setcolors,/sys
thisimg=FLUXBALUNSIMAGE[*,30:150]
thisimg=thisimg/max(thisimg)*242.
thisimg[where(fluximage[*,30:150] eq 0)]=253
thisimg[0,0]=254.
plot_image,thisimg,scale=[xbin,1.],orig=[0.,-60.],ytit='Latitude',xtit='Years Since '+anytim(mintim*3600.*24.*365.,/vms,/date),/nosq
sharpcorners
;plot_image,FLUXBALUNSIMAGE[*,30:150],/nosq,scale=[xbin,1.],orig=[0.,-60.]
;oplot,tim1pd[wgt22],arstr_arr[wgt22].hglat,ps=8,color=!red
;oplot,tim1pd[wwlsg[wle60]],arstr_arr[wwlsg[wle60]].hglat,ps=8,color=!green
loadct,1
;color_table,[min(FLUXBALUNSIMAGE),max(FLUXBALUNSIMAGE)],tit='<'+textoidl('\Phi')+' Imb.>'
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_fluximb.eps '+paper2path+'butterfly_2wk_fluximb.pdf';spawn,'ps2pdf -dEPSCrop '+paper2path+'butterfly_2wk_fluximb.eps '+paper2path+'butterfly_2wk_fluximb.pdf'

;unsigned flux balance in all regions
setplotenv,file=paper2path+'butterfly_2wk_fluximb_net.eps',/ps
!p.multi=[0,1,2]
loadct,0
thisimg0=alog10((fluximage[*,30:149]*1d16 < 2d23) + 1.)*(-1.) & thisimg0[where(thisimg0 eq 0)]=-21.8348;mean(minmax(thisimg0 < max(thisimg0)));(max(thisimg0)-(min(thisimg0)>0.))/2.;6000000./2.
plot_image,thisimg0,scale=[xbin,1.],orig=[0.,-60.],ytit='Latitude',/nosq,ymarg=[0,1],xtickname=strarr(10)+' '
hline,0,color=0
thisimg=(abs(fluxposimage-fluxnegimage)/(fluximage))[*,30:149]*(1.) & thisimg[where(finite(thisimg) ne 1)]=.5
plot_image,thisimg,scale=[xbin,1.],orig=[0.,-60.],ytit='Latitude',xtit='Years Since 1997',/nosq,ymarg=[3,0]
hline,0,color=0
color_table, minmax((fluximage[*,30:149] < 6000000)*1d16),[.7,.95],[.9,.95],rr=rr0,gg=gg0,bb=bb0,title=textoidl('\Phi_{tot}')+' [Mx]'
color_table, minmax(thisimg),[.7,.95],[.42,.47],rr=rr0,gg=gg0,bb=bb0,title=textoidl('|\Phi_{+}+\Phi_{-}|/(\Phi_{tot})')+' [Mx]'
sharpcorners
;setcolors,/sys
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_fluximb_net.eps '+paper2path+'butterfly_2wk_fluximb_net.pdf'


stop

;PLOT hglat flares complexity----------------------------------------->
plotsym,0,.5,/fill
setcolors,/sys

plot,tim1pd,arstr_arr.hglat,ps=8,yran=[-60,60],xran=[0,11.2],/xsty,/ysty,xtit='years since 19970101',ytit='hg latitude',tit='R=10^22mx G=15k WLSG B=MX Flares'
oplot,tim1pd,arstr_arr.hglat,ps=8,color=150

plotsym,0,.7,/fill
wgt22=where(arstr_arr.bflux*1d16 gt 1d22)
oplot,tim1pd[wgt22],arstr_arr[wgt22].hglat,ps=8,color=!red

wwlsg=where(arstr_arr.nlstr.wlsg gt 15000)
oplot,tim1pd[wwlsg],arstr_arr[wwlsg].hglat,ps=8,color=!green
wle60=where(arstr_arr[wwlsg].hglon le 60.)
oplot,tim1pd[wwlsg[wle60]],arstr_arr[wwlsg[wle60]].hglat,ps=8,color=!forest

plotsym,0,.7
wmx=where(strmid(ASTRPOS.GOES_CLASS,0,1) eq 'M' or strmid(ASTRPOS.GOES_CLASS,0,1) eq 'X')
window_capture,file=datapath+'hglat_wlsg_mxflares_plot0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,(FLRPOSSTR[wmx].hgpos)[1,*],ps=8,color=0
window_capture,file=datapath+'hglat_wlsg_mxflares_plot1',/png

stop

;IMAGE hglat flares complexity----------------------------------------->
loadct,0
plot_image,(rvalimage)^.2,/nosq,origin=[0,-90],scale=[xbin,ybin],xtit='years since 19970101',ytit='HG Latitude',tit='R value' 
plotsym,0,.5,/fill
setcolors,/sys
window_capture,file=datapath+'hglat_rvalue_mxflares_image0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,(FLRPOSSTR[wmx].hgpos)[1,*],ps=8,color=!red
window_capture,file=datapath+'hglat_rvalue_mxflares_image1',/png

loadct,0
plot_image,(wlsgimage)^.2,/nosq,origin=[0,-90],scale=[xbin,ybin],xtit='years since 19970101',ytit='HG Latitude',tit='WLSG' 
plotsym,0,.5,/fill
setcolors,/sys
window_capture,file=datapath+'hglat_wlsg_mxflares_image0',/png

oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,(FLRPOSSTR[wmx].hgpos)[1,*],ps=8,color=!red
window_capture,file=datapath+'hglat_wlsg_mxflares_image1',/png
stop

;PLOT flux oscillation in N/S Hemispheres------------------------------->
;!p.multi=[0,1,2]
!p.multi=[2,1,2]
plot,total(flrmximage,2),xsty=4,ysty=4,color=!blue
oplot,total(flrindeximage,2)/max(total(flrindeximage,2))*max(total(flrmximage,2)),color=!cyan
!p.multi=[2,1,2]
plot,total(fluximage,2),ytit='bflux G/Mm^2'
oplot,total(fluximage[*,0:89],2),color=!green
oplot,total(fluximage[*,90:179],2),color=!red

yarr=findgen((size(fluximage))[2])
bfluxlatcent=fltarr((size(fluximage))[1])
for i=0,n_elements(bfluxlatcent)-1 do bfluxlatcent[i]=total(fluximage[i,*]*yarr)/total(fluximage[i,*])
yarrs=yarr[0:89]
bfluxlatcents=bfluxlatcent
for i=0,n_elements(bfluxlatcents)-1 do bfluxlatcents[i]=total(fluximage[i,0:89]*yarrs)/total(fluximage[i,0:89])
yarrn=yarr[90:179]
bfluxlatcentn=bfluxlatcent
for i=0,n_elements(bfluxlatcentn)-1 do bfluxlatcentn[i]=total(fluximage[i,90:179]*yarrn)/total(fluximage[i,90:179])
bfluxlatcent[where(finite(bfluxlatcent) ne 1)]=0 & bfluxlatcents[where(finite(bfluxlatcents) ne 1)]=0 & bfluxlatcentn[where(finite(bfluxlatcentn) ne 1)]=0

!p.multi=[1,1,2]
plot,bfluxlatcent,ytit='hg lat'
oplot,bfluxlatcents,color=!green
oplot,bfluxlatcentn,color=!red
window_capture,file=datapath+'oscillation_2wkbin_bflux_latitude_all',/png

;PLOT flux oscillation in N/S Hemispheres SMOOTHED---------------------->
erase
!p.multi=[2,1,2]
plot,smooth(total(flrmximage,2),4),xsty=4,ysty=4,color=!blue
oplot,smooth(total(flrindeximage,2),4)/max(smooth(total(flrindeximage,2),4))*max(smooth(total(flrmximage,2),4)),color=!cyan

!p.multi=[2,1,2]
plot,smooth(total(fluximage,2),4),ytit='bflux G Mm^2',/noerase
oplot,smooth(total(fluximage[*,0:89],2),4),color=!green
oplot,smooth(total(fluximage[*,90:179],2),4),color=!red
!p.multi=[1,1,2]
plot,(bfluxlatcent),ytit='hg lat'
oplot,(bfluxlatcents),color=!green
oplot,(bfluxlatcentn),color=!red
window_capture,file=datapath+'oscillation_2wkbin_smth_bflux_latitude_all',/png

;PLOT flux oscillation in N/S Hemispheres from bipolar flux image------->
yarr=findgen((size(fluxbipolimage))[2])
bfluxlatcent=fltarr((size(fluxbipolimage))[1])
for i=0,n_elements(bfluxlatcent)-1 do bfluxlatcent[i]=total(fluxbipolimage[i,*]*yarr)/total(fluxbipolimage[i,*])
yarrs=yarr[0:89]
bfluxlatcents=bfluxlatcent
for i=0,n_elements(bfluxlatcents)-1 do bfluxlatcents[i]=total(fluxbipolimage[i,0:89]*yarrs)/total(fluxbipolimage[i,0:89])
yarrn=yarr[90:179]
bfluxlatcentn=bfluxlatcent
for i=0,n_elements(bfluxlatcentn)-1 do bfluxlatcentn[i]=total(fluxbipolimage[i,90:179]*yarrn)/total(fluxbipolimage[i,90:179])
bfluxlatcent[where(finite(bfluxlatcent) ne 1)]=0 & bfluxlatcents[where(finite(bfluxlatcents) ne 1)]=0 & bfluxlatcentn[where(finite(bfluxlatcentn) ne 1)]=0

plot,smooth(total(fluxbipolimage,2),4),ytit='bflux G/Mm^2'
oplot,smooth(total(fluxbipolimage[*,0:89],2),4),color=!green
oplot,smooth(total(fluxbipolimage[*,90:179],2),4),color=!red
plot,(bfluxlatcent),ytit='hg lat'
oplot,(bfluxlatcents),color=!green
oplot,(bfluxlatcentn),color=!red
window_capture,file=datapath+'oscillation_2wkbin_smth_bflux_latitude_bipolar',/png

stop

;Compare bipole and unipole flux over SC.
setplotenv,file=paper2path+'bipole_unipole_flux_compare.eps',/ps,xs=18,ys=14
plot,findgen(292)*xbin,smooth(total(fluxbipolimage[*,90:179],2),5),ytit='North Hemisphere Multipolar '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']',position=[.1,.5,.9,.95],xtickname=strarr(10)+' ',/xsty
;oplot,deriv(findgen(292)*xbin,smooth(total(fluxbipolimage[*,90:179],2),11)),color=150
plot,findgen(292)*xbin,smooth(total(fluxunipolimage[*,90:179],2),5),lines=2,/xsty,position=[.1,.5,.9,.95],yticklen=.001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' ',/ysty,yran=[0,.5d7],/noerase
axis,yaxis=1,yran=[0,.5d7],ytit='Unipolar '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']'
legend,['Multipolar '+textoidl('\Phi'),'Unipolar '+textoidl('\Phi')],lines=[0,2],/right,/top
;legend,['Multipolar '+textoidl('\Phi'),'Multipolar '+textoidl('d\Phi/dt'),'Unipolar '+textoidl('\Phi')],lines=[0,0,2],color=[0,150,0]

plot,findgen(292)*xbin,smooth(total(fluxbipolimage[*,0:89],2),5),ytit='South Hemisphere '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']',position=[.1,.1,.9,.5],/noerase,/xsty,xtit='Years Since 3-Feb-2007'
;oplot,deriv(findgen(292)*xbin,smooth(total(fluxbipolimage[*,0:89],2),11)),color=150
plot,findgen(292)*xbin,smooth(total(fluxunipolimage[*,0:89],2),5),lines=2,position=[.1,.1,.9,.5],/ysty,yran=[0,.6d7],/noerase,/xsty,yticklen=.001,ytickname=strarr(10)+' '
axis,yaxis=1,yran=[0,.6d7],ytit='Unipolar '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']'
;legend,['Multipolar '+textoidl('\Phi'),'Multipolar '+textoidl('d\Phi/dt'),'Unipolar '+textoidl('\Phi')],lines=[0,0,2],color=[0,150,0]
closeplotenv & spawn,'convert '+paper2path+'bipole_unipole_flux_compare.eps '+paper2path+'bipole_unipole_flux_compare.pdf'

ccn=cross_corr2(rebin(smooth(total(fluxunipolimage[*,0:89],2),7),292,2),rebin(smooth(total(fluxbipolimage[*,0:89],2),7),292,2),[20,1],cx,cy,miss=wmiss,/report)
ccs=cross_corr2(rebin(smooth(total(fluxunipolimage[*,90:179],2),7),292,2),rebin(smooth(total(fluxbipolimage[*,90:179],2),7),292,2),[20,1],cx,cy,miss=wmiss,/report)

;Compare total hemispheric Flux to the .
;setplotenv,file=paper2path+'bipole_unipole_flux_compare.eps',/ps,xs=18,ys=14
;plot,findgen(292)*xbin,smooth(total(fluxbipolimage[*,90:179],2),5),ytit='North Hemisphere Multipolar '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']',position=[.1,.5,.9,.95],xtickname=strarr(10)+' ',/xsty
;;oplot,deriv(findgen(292)*xbin,smooth(total(fluxbipolimage[*,90:179],2),11)),color=150
;plot,findgen(292)*xbin,smooth(total(fluxunipolimage[*,90:179],2),5),lines=2,/xsty,position=[.1,.5,.9,.95],yticklen=.001,xtickname=strarr(10)+' ',ytickname=strarr(10)+' ',/ysty,yran=[0,.5d7],/noerase
;axis,yaxis=1,yran=[0,.5d7],ytit='Unipolar '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']'
;legend,['Multipolar '+textoidl('\Phi'),'Unipolar '+textoidl('\Phi')],lines=[0,2],/right,/top
;;legend,['Multipolar '+textoidl('\Phi'),'Multipolar '+textoidl('d\Phi/dt'),'Unipolar '+textoidl('\Phi')],lines=[0,0,2],color=[0,150,0]

;plot,findgen(292)*xbin,smooth(total(fluxbipolimage[*,0:89],2),5),ytit='South Hemisphere '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']',position=[.1,.1,.9,.5],/noerase,/xsty,xtit='Years Since 3-Feb-2007'
;;oplot,deriv(findgen(292)*xbin,smooth(total(fluxbipolimage[*,0:89],2),11)),color=150
;plot,findgen(292)*xbin,smooth(total(fluxunipolimage[*,0:89],2),5),lines=2,position=[.1,.1,.9,.5],/ysty,yran=[0,.6d7],/noerase,/xsty,yticklen=.001,ytickname=strarr(10)+' '
;axis,yaxis=1,yran=[0,.6d7],ytit='Unipolar '+textoidl('\Phi')+' ['+textoidl('Mm G^2')+']'
;;legend,['Multipolar '+textoidl('\Phi'),'Multipolar '+textoidl('d\Phi/dt'),'Unipolar '+textoidl('\Phi')],lines=[0,0,2],color=[0,150,0]
;closeplotenv & spawn,'convert '+paper2path+'bipole_unipole_flux_compare.eps '+paper2path+'bipole_unipole_flux_compare.pdf'


;PLOT flux oscillation in N/S Hemispheres from GT 10^22 Mx image-------->

yarr=findgen((size(fluxgt22image))[2])
bfluxlatcent=fltarr((size(fluxgt22image))[1])
for i=0,n_elements(bfluxlatcent)-1 do bfluxlatcent[i]=total(fluxgt22image[i,*]*yarr)/total(fluxgt22image[i,*])
yarrs=yarr[0:89]
bfluxlatcents=bfluxlatcent
for i=0,n_elements(bfluxlatcents)-1 do bfluxlatcents[i]=total(fluxgt22image[i,0:89]*yarrs)/total(fluxgt22image[i,0:89])
yarrn=yarr[90:179]
bfluxlatcentn=bfluxlatcent
for i=0,n_elements(bfluxlatcentn)-1 do bfluxlatcentn[i]=total(fluxgt22image[i,90:179]*yarrn)/total(fluxgt22image[i,90:179])
bfluxlatcent[where(finite(bfluxlatcent) ne 1)]=0 & bfluxlatcents[where(finite(bfluxlatcents) ne 1)]=0 & bfluxlatcentn[where(finite(bfluxlatcentn) ne 1)]=0

plot,smooth(total(fluxgt22image,2),4),ytit='bflux G/Mm^2'
oplot,smooth(total(fluxgt22image[*,0:89],2),4),color=!green
oplot,smooth(total(fluxgt22image[*,90:179],2),4),color=!red
plot,(bfluxlatcent),ytit='hg lat'
oplot,(bfluxlatcents),color=!green
oplot,(bfluxlatcentn),color=!red
window_capture,file=datapath+'oscillation_2wkbin_smth_bflux_latitude_gt22',/png

stop

;PLOT flux oscillation in N/S Hemispheres from GT 10^23 Mx image-------->

yarr=findgen((size(fluxgt23image))[2])
bfluxlatcent=fltarr((size(fluxgt23image))[1])
for i=0,n_elements(bfluxlatcent)-1 do bfluxlatcent[i]=total(fluxgt23image[i,*]*yarr)/total(fluxgt23image[i,*])
yarrs=yarr[0:89]
bfluxlatcents=bfluxlatcent
for i=0,n_elements(bfluxlatcents)-1 do bfluxlatcents[i]=total(fluxgt23image[i,0:89]*yarrs)/total(fluxgt23image[i,0:89])
yarrn=yarr[90:179]
bfluxlatcentn=bfluxlatcent
for i=0,n_elements(bfluxlatcentn)-1 do bfluxlatcentn[i]=total(fluxgt23image[i,90:179]*yarrn)/total(fluxgt23image[i,90:179])
bfluxlatcent[where(finite(bfluxlatcent) ne 1)]=0 & bfluxlatcents[where(finite(bfluxlatcents) ne 1)]=0 & bfluxlatcentn[where(finite(bfluxlatcentn) ne 1)]=0

plot,smooth(total(fluxgt23image,2),4),ytit='bflux G/Mm^2'
oplot,smooth(total(fluxgt23image[*,0:89],2),4),color=!green
oplot,smooth(total(fluxgt23image[*,90:179],2),4),color=!red
plot,(bfluxlatcent),ytit='hg lat'
oplot,(bfluxlatcents),color=!green
oplot,(bfluxlatcentn),color=!red
window_capture,file=datapath+'oscillation_2wkbin_smth_bflux_latitude_gt23',/png

;PLOT flux, complexity, flare numbers, north and south hemispheres------------------>
!p.multi=[0,1,3]
plot,smooth(total(flrmximage,2),4),ytit='All Bk=FLR R=Bflx Bl=Rval G=Wlsg',ymargin=[0,0],chars=3, pos=[.1,.05,.9,.35]
oplot,smooth(total(flrindeximage,2),4)/3.,color=150
oplot,smooth(total(fluximage,2),4)/3d7,color=!red
oplot,smooth(total(rvalimage,2),4)/1d13,color=!blue
oplot,smooth(total(wlsgimage,2),4)/5d4,color=!green
plot,smooth(total(flrmximage[*,0:89],2),4),ytit='SOUTH',ymargin=[0,0],chars=3, pos=[.1,.35,.9,.65]
oplot,smooth(total(flrindeximage[*,0:89],2),4)/3.,color=150
oplot,smooth(total(fluximage[*,0:89],2),4)/3d7,color=!red
oplot,smooth(total(rvalimage[*,0:89],2),4)/1d13,color=!blue
oplot,smooth(total(wlsgimage[*,0:89],2),4)/5d4,color=!green
plot,smooth(total(flrmximage[*,90:179],2),4),ytit='NORTH',ymargin=[0,0],chars=3, pos=[.1,.65,.9,.95]
oplot,smooth(total(flrindeximage[*,90:179],2),4)/3.,color=150
oplot,smooth(total(fluximage[*,90:179],2),4)/3d7,color=!red
oplot,smooth(total(rvalimage[*,90:179],2),4)/1d13,color=!blue
oplot,smooth(total(wlsgimage[*,90:179],2),4)/5d4,color=!green
window_capture,file=datapath+'oscillation_2wkbin_smth_magprop_vs_flr_all_n_s',/png

stop

;unsigned flux in unipolar regions-------------------------------------->
loadct,1
plot_image,FLUXUNIPOLIMAGE,/nosqr

;unsigned flux balance in unipolar regions
loadct,0
plot_image,FLUXBALUNIPOLIMAGE,/nosqr

;flux emergence and decay regions
loadct,0
plot_image,BFLUXEMERGEIMAGE < max(BFLUXEMERGEIMAGE)*.001 > (-max(BFLUXEMERGEIMAGE)*.001),/nosqr

;area 
loadct,0
plot_image,(areaimage[*,30:150])^.5,/nosq,orig=[0,-60],scale=[xbin,ybin]
setcolors,/sys
oplot,FLRPOSSTR[wmx].tim/(3600.*24.*365)-mintim,(FLRPOSSTR[wmx].hgpos)[1,*],ps=8,color=!red

;2d colar image!!! wlsg over flux SMOOTHED!!
setplotenv,file=paper2path+'butterfly_2wk_2dcolor_flux_wlsg_smth.eps',/ps,xs=15,ys=12
!p.multi=0
xx=fltarr(292,180,3)
xwlsgimage=smooth(wlsgimage,[3,3])
xfluximage=smooth(fluximage,[3,3])
xx[*,*,1]=xwlsgimage/max(xwlsgimage)                        
xx[*,*,0]=xfluximage/max(xfluximage)
xxblankvals=fltarr(292,180)
wblanks=where(xwlsgimage eq 0 and xfluximage eq 0)
xxblankvals[wblanks]=1 & xx[*,*,2]=xxblankvals
xx0=xx[*,*,0] & xx0[wblanks]=1 & xx[*,*,0]=xx0
xx1=xx[*,*,1] & xx1[wblanks]=1 & xx[*,*,1]=xx1
plot_image,(abs(xx)+.001)^.3,true=3,position=[.1,.1,.9,.9],origin=[0,-90],scale=[xbin,1.],yrange=[-60,60],ytit='Latitude',xtit='Years since 3-Feb-1997';,/NOADJ
ct2d=fltarr(255,255,3)
ctg=rebin(findgen(255),255,255)
ctr=rot(ctg,-90)
ct2d[*,*,0]=ctr
ct2d[*,*,1]=ctg
plot_image,(abs(ct2d)+.001)^.3,true=3,position=[.70,.7,.89,.89],/noerase,ytit='Flux',xtit='WLsg',scale=[max(xwlsgimage)/254.,max(xfluximage)/254.],xticks=1,xtickv=[2d5,4d5],yticks=1,ytickv=[1d8,2d8],xticklen=1,yticklen=1,xminor=0,yminor=0,chars=2
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_2dcolor_flux_wlsg_smth.eps '+paper2path+'butterfly_2wk_2dcolor_flux_wlsg_smth.pdf'

;2d colar image!!! wlsg over flux
setplotenv,file=paper2path+'butterfly_2wk_2dcolor_flux_wlsg.eps',/ps,xs=15,ys=12
!p.multi=0
xx=fltarr(292,180,3)
xwlsgimage=wlsgimage
xfluximage=fluximage
xx[*,*,1]=xwlsgimage/max(xwlsgimage)                        
xx[*,*,0]=xfluximage/max(xfluximage)
xxblankvals=fltarr(292,180)
wblanks=where(xwlsgimage eq 0 and xfluximage eq 0)
xxblankvals[wblanks]=1 & xx[*,*,2]=xxblankvals
xx0=xx[*,*,0] & xx0[wblanks]=1 & xx[*,*,0]=xx0
xx1=xx[*,*,1] & xx1[wblanks]=1 & xx[*,*,1]=xx1
plot_image,(abs(xx)+.001)^.3,true=3,position=[.1,.1,.9,.9],origin=[0,-90],scale=[xbin,1.],yrange=[-60,60],ytit='Latitude',xtit='Years since 3-Feb-1997';,/NOADJ
ct2d=fltarr(255,255,3)
ctg=rebin(findgen(255),255,255)
ctr=rot(ctg,-90)
ct2d[*,*,0]=ctr
ct2d[*,*,1]=ctg
plot_image,(abs(ct2d)+.001)^.3,true=3,position=[.70,.7,.89,.89],/noerase,ytit='Flux',xtit='WLsg',scale=[max(xwlsgimage)/254.,max(xfluximage)/254.],xticks=1,xtickv=[2d5,4d5],yticks=1,ytickv=[1d8,2d8],xticklen=1,yticklen=1,xminor=0,yminor=0,chars=2
closeplotenv & spawn,'convert '+paper2path+'butterfly_2wk_2dcolor_flux_wlsg.eps '+paper2path+'butterfly_2wk_2dcolor_flux_wlsg.pdf'


stop


;-------------------------------------------------------------------------------------------------->
;REBIN TO 2months---------------------------------------------------------------------------------->

if not keyword_set(res1tore) then begin
;Run latitude images
ybin=1. & xbin=0.166667
xx=smart_bin2d(xbin=xbin,ybin=ybin,xran=xran,yran=yran,/years,struct=arstr_arr, $
	flrstruct=astrpos,flrtimarr=FLRPOSSTR.tim/3600./24./365.-mintim, flrposstr=FLRPOSSTR, $
	/latitude,timarr=tim1pd,path=datapath,outfile=outfile,filemod='_2month_bin')
restore,outfile,/ver
endif else restore,datapath+'butterfly_diags_lat_2month_bin.sav'

stop

;2d colar image!!! wlsg over flux
setplotenv,file=paper2path+'butterfly_2dcolor_flux_wlsg.eps',/ps,xs=15,ys=12
!p.multi=0
xx=fltarr(67,180,3)
xx[*,*,1]=wlsgimage/max(wlsgimage)                        
xx[*,*,0]=fluximage/max(fluximage)
xxblankvals=fltarr(67,180)
wblanks=where(wlsgimage eq 0 and fluximage eq 0)
xxblankvals[wblanks]=1 & xx[*,*,2]=xxblankvals
xx0=xx[*,*,0] & xx0[wblanks]=1 & xx[*,*,0]=xx0
xx1=xx[*,*,1] & xx1[wblanks]=1 & xx[*,*,1]=xx1
plot_image,(abs(xx)+.001)^.3,true=3,position=[.1,.1,.9,.9],origin=[0,-90],scale=[xbin,1.],yrange=[-60,60],ytit='Latitude',xtit='Years since 3-Feb-1997';,/NOADJ
ct2d=fltarr(255,255,3)
ctg=rebin(findgen(255),255,255)
ctr=rot(ctg,-90)
ct2d[*,*,0]=ctr
ct2d[*,*,1]=ctg
plot_image,fltarr(255,255,3),thick=1,true=3,position=[.70,.7,.89,.89],/noerase,ytit='Flux',xtit='WLsg',scale=[max(wlsgimage)/254.,max(fluximage)/254.],xticks=1,xtickv=[2d5,4d5],yticks=1,ytickv=[1d8,2d8],xticklen=1,yticklen=1,xminor=0,yminor=0,chars=2
plot_image,(abs(ct2d)+.001)^.3,thick=1,true=3,position=[.70,.7,.89,.89],/noerase,ytit='Flux',xtit='WLsg',scale=[max(wlsgimage)/254.,max(fluximage)/254.],xticks=1,xtickv=[2d5,4d5],yticks=1,ytickv=[1d8,2d8],xticklen=1,yticklen=1,xminor=0,yminor=0,chars=2
closeplotenv & spawn,'convert '+paper2path+'butterfly_2dcolor_flux_wlsg.eps '+paper2path+'butterfly_2dcolor_flux_wlsg.pdf'

;2d colar image!!! wlsg over flux SMOOTHED!!
setplotenv,file=paper2path+'butterfly_2dcolor_flux_wlsg_smth.eps',/ps,xs=15,ys=12
!p.multi=0
xx=fltarr(67,180,3)
xwlsgimage=smooth(wlsgimage,[3,3])
xfluximage=smooth(fluximage,[3,3])
xx[*,*,1]=xwlsgimage/max(xwlsgimage)                        
xx[*,*,0]=xfluximage/max(xfluximage)
xxblankvals=fltarr(67,180)
wblanks=where(wlsgimage eq 0 and fluximage eq 0)
xxblankvals[wblanks]=1 & xx[*,*,2]=xxblankvals
xx0=xx[*,*,0] & xx0[wblanks]=1 & xx[*,*,0]=xx0
xx1=xx[*,*,1] & xx1[wblanks]=1 & xx[*,*,1]=xx1
plot_image,(abs(xx)+.001)^.3,true=3,position=[.1,.1,.9,.9],origin=[0,-90],scale=[xbin,1.],yrange=[-60,60],ytit='Latitude',xtit='Years since 3-Feb-1997';,/NOADJ
ct2d=fltarr(255,255,3)
ctg=rebin(findgen(255),255,255)
ctr=rot(ctg,-90)
ct2d[*,*,0]=ctr
ct2d[*,*,1]=ctg
plot_image,(abs(ct2d)+.001)^.3,true=3,position=[.70,.7,.89,.89],/noerase,ytit='Flux',xtit='WLsg',scale=[max(xwlsgimage)/254.,max(xfluximage)/254.],xticks=1,xtickv=[2d5,4d5],yticks=1,ytickv=[1d8,2d8],xticklen=1,yticklen=1,xminor=0,yminor=0,chars=2
closeplotenv & spawn,'convert '+paper2path+'butterfly_2dcolor_flux_wlsg_smth.eps '+paper2path+'butterfly_2dcolor_flux_wlsg_smth.pdf'


;PLOT flux oscillation in N/S Hemispheres------------------------------->
!p.multi=[0,1,2]
plot,total(fluximage,2),ytit='bflux G/Mm^2'
oplot,total(fluximage[*,0:89],2),color=!green
oplot,total(fluximage[*,90:179],2),color=!red

yarr=findgen((size(fluximage))[2])
bfluxlatcent=fltarr((size(fluximage))[1])
for i=0,n_elements(bfluxlatcent)-1 do bfluxlatcent[i]=total(fluximage[i,*]*yarr)/total(fluximage[i,*])
yarrs=yarr[0:89]
bfluxlatcents=bfluxlatcent
for i=0,n_elements(bfluxlatcents)-1 do bfluxlatcents[i]=total(fluximage[i,0:89]*yarrs)/total(fluximage[i,0:89])
yarrn=yarr[90:179]
bfluxlatcentn=bfluxlatcent
for i=0,n_elements(bfluxlatcentn)-1 do bfluxlatcentn[i]=total(fluximage[i,90:179]*yarrn)/total(fluximage[i,90:179])

plot,bfluxlatcent,ytit='hg lat'
oplot,bfluxlatcents,color=!green
oplot,bfluxlatcentn,color=!red





















end