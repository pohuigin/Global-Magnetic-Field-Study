function smart_bin2d, xarr, yarr, xbin=xbin, ybin=ybin, xrange=xrange, yrange=yrange, years=years, $
	struct=instruct, flrstruct=inflrstruct, timarr=timarr, flrtimarr=flrtimarr, flrposstr=flrposstr, latitude=latitude, longitude=longitude, longarr=longarr, $
	path=path, outfile=outfile,filemod=filemod

if n_elements(path) lt 1 then datapath='~/science/data/cycle23_sav/smart_paper_2/' else datapath=path

;if n_elements(instruct) lt 1 then struct
struct=instruct
if n_elements(inflrstruct) gt 0 then flrstruct=inflrstruct

wlt60=where(struct.hglon lt 60 and struct.hglon gt -60)
struct=struct[wlt60]
timarr=timarr[wlt60]

;tarr=xarr
;darr=yarr
;bin2=2
;bin1=.083
;min1=0
;min2=0
;xx=hist_2d(tarr,darr,min1=min1,max1=11.2,bin1=bin1,min2=min2,max2=360-bin2,bin2=bin2)
;xs=findgen((size(xx))[1])
;xs=xs*bin1+min1
;ys=findgen((size(xx))[2])
;ys=ys*bin2+min2

nx=round(abs(xrange[0]-xrange[1])/float(xbin))
ny=round(abs(yrange[0]-yrange[1])/float(ybin))

;window,xs=nx,ys=ny

fluximage=fltarr(nx,ny)
fluxgt22image=fltarr(nx,ny)
fluxgt23image=fltarr(nx,ny)
fluxbipolimage=fltarr(nx,ny)
areaimage=fltarr(nx,ny)
fluxbalimage=fltarr(nx,ny)
fluxbalunsimage=fltarr(nx,ny)
fluxunipolimage=fltarr(nx,ny)
fluxbalunipolimage=fltarr(nx,ny)
fluxbalbipolimage=fltarr(nx,ny)
rvalimage=fltarr(nx,ny)
wlsgimage=fltarr(nx,ny)
rvalstarimage=fltarr(nx,ny)
wlsgstarimage=fltarr(nx,ny)
bfluxemergeimage=fltarr(nx,ny)

;PAPER 2: FLUXBALANCE (new plots)
fluxsignedimage=fltarr(nx,ny)
fluxposimage=fltarr(nx,ny)
fluxnegimage=fltarr(nx,ny)
fluxfracimage=fltarr(nx,ny)
fluxfracavgimage=fltarr(nx,ny)
fluxunipolpos=fltarr(nx,ny)
fluxunipolneg=fltarr(nx,ny)
fluxbipolpos=fltarr(nx,ny)
fluxbipolneg=fltarr(nx,ny)

;PAPER 2: dFLUX/dt (new plots)
fluxuniemerge=fltarr(nx,ny)
fluxuniemergeavg=fltarr(nx,ny)

;Flare stuff
flrindeximage=fltarr(nx,ny)
flrnumimage=fltarr(nx,ny)
flrmximage=fltarr(nx,ny)

if keyword_set(years) then secpyr=3600.*24.*365. else secpyr=1.

if n_elements(timarr) eq 0 then tim=anytim(struct.time)/secpyr else tim=timarr ;-min(anytim(struct.time)/3600./24./365.)
if keyword_set(latitude) then lat=struct.hglat else begin
	if n_elements(longarr) gt 0 then lat=longarr else lat=struct.carlon
endelse

if n_elements(inflrstruct) gt 0 then begin
	if n_elements(flrtimarr) eq 0 then flrtim=anytim(flrstruct.time)/secpyr else flrtim=flrtimarr
	if keyword_set(latitude) then flrlat=(flrposstr.hgpos)[1,*] else flrlat=(flrposstr.carlon)
endif

xarr=findgen(nx)/(nx-1.)*abs(xrange[0]-xrange[1]+xbin)+xrange[0]
yarr=findgen(ny)/(ny-1.)*abs(yrange[0]-yrange[1]+ybin)+yrange[0]
for i=0,nx-1 do begin

;if i ge 214 then stop
	
	wthistim=where(tim ge (i*xbin+xrange[0]) and tim lt ((i+1.)*xbin+xrange[0]))
	if n_elements(inflrstruct) gt 0 then wflrtim=where(flrtim ge (i*xbin+xrange[0]) and flrtim lt ((i+1.)*xbin+xrange[0]))
	
	for j=0,ny-1 do begin
	if wthistim[0] ne -1 then begin
			wthisy=where(lat[wthistim] ge (j*ybin+yrange[0]) and lat[wthistim] lt ((j+1.)*ybin+yrange[0]))
			if wthisy[0] ne -1 then begin
		
				thisstruct=struct[wthistim[wthisy]]
		
				fluximage[i,j]=total(thisstruct.bflux)
				areaimage[i,j]=total(thisstruct.area)
				fluxbalimage[i,j]=mean((thisstruct.bfluxpos-thisstruct.bfluxneg)/thisstruct.bflux)
				fluxbalunsimage[i,j]=mean(abs(thisstruct.bfluxpos-thisstruct.bfluxneg)/thisstruct.bflux)
				wuni=where(strlowcase((thisstruct.type)[0]) eq 'u')
				if wuni[0] ne -1 then fluxunipolimage[i,j]=total(thisstruct[wuni].bflux)
				if wuni[0] ne -1 then fluxbalunipolimage[i,j]=mean((thisstruct[wuni].bfluxpos-thisstruct[wuni].bfluxneg)/thisstruct[wuni].bflux)
				rvalimage[i,j]=total(thisstruct.NLSTR.rval)
				wlsgimage[i,j]=total(thisstruct.NLSTR.wlsg)
				rvalstarimage[i,j]=total(thisstruct.NLSTR.r_star)
				wlsgstarimage[i,j]=total(thisstruct.NLSTR.wlsg_star)
				bfluxemergeimage[i,j]=total(thisstruct.BFLUXEMRG)
				if wuni[0] ne -1 then fluxuniemerge[i,j]=total(thisstruct[wuni].BFLUXEMRG)
				if wuni[0] ne -1 then begin
					if (where(thisstruct[wuni].BFLUXEMRG gt 0))[0] ne -1 then nuniemrg=n_elements(where(thisstruct[wuni].BFLUXEMRG gt 0)) else nuniemrg=0.
					if (where(thisstruct[wuni].BFLUXEMRG lt 0))[0] ne -1 then nunidecy=n_elements(where(thisstruct[wuni].BFLUXEMRG lt 0)) else nunidecy=0.
					if (nuniemrg+nunidecy) eq 0 then fluxuniemergeavg[i,j]=0. else fluxuniemergeavg[i,j]=(nuniemrg-nunidecy)/float(nuniemrg+nunidecy)
				endif
				wbi=where(strlowcase((thisstruct.type)[0]) eq 'm')
				if wbi[0] ne -1 then fluxbipolimage[i,j]=total(thisstruct[wbi].bflux)
				if wbi[0] ne -1 then fluxbalbipolimage[i,j]=mean((thisstruct[wbi].bfluxpos-thisstruct[wbi].bfluxneg)/thisstruct[wbi].bflux)
				wgt22=where(thisstruct.bflux*1d16 gt 1d22)
				if wgt22[0] ne -1 then fluxgt22image[i,j]=total(thisstruct[wgt22].bflux)
				wgt23=where(thisstruct.bflux*1d16 gt 1d23)
				if wgt23[0] ne -1 then fluxgt23image[i,j]=total(thisstruct[wgt23].bflux)

				fluxsignedimage[i,j]=total(thisstruct.bfluxpos)-abs(total(thisstruct.bfluxneg))
				fluxposimage[i,j]=total(thisstruct.bfluxpos)
				fluxnegimage[i,j]=abs(total(thisstruct.bfluxneg))
				fluxfracimage[i,j]=(total(thisstruct.bfluxpos) < abs(total(thisstruct.bfluxneg)))/(total(thisstruct.bfluxpos) > abs(total(thisstruct.bfluxneg)))
				fluxfracavgimage[i,j]=mean((thisstruct.bfluxpos < abs(thisstruct.bfluxneg))/(thisstruct.bfluxpos > abs(thisstruct.bfluxneg)))
				if wuni[0] ne -1 then fluxunipolpos[i,j]=total(thisstruct[wuni].bfluxpos)
				if wuni[0] ne -1 then fluxunipolneg[i,j]=total(thisstruct[wuni].bfluxneg)
				
				if wbi[0] ne -1 then fluxbipolpos[i,j]=total(thisstruct[wbi].bfluxpos)
				if wbi[0] ne -1 then fluxbipolneg[i,j]=total(thisstruct[wbi].bfluxneg)

			endif
		endif
;stop
		if n_elements(inflrstruct) gt 0 then begin
			if wflrtim[0] ne -1 then begin
				wflry=where(flrlat[wflrtim] ge (j*ybin+yrange[0]) and flrlat[wflrtim] lt ((j+1.)*ybin+yrange[0]))
				if wflry[0] ne -1 then begin
					thisflrstr=flrstruct[wflrtim[wflry]]
					if keyword_set(years) then tau=xbin*365. else tau=xbin/3600./24.
					
					flrindeximage[i,j]=smart_goes_findex(thisflrstr, tau)
					flrnumimage[i,j]=n_elements(thisflrstr)
					flrmximage[i,j]=n_elements(where(strmid(thisflrstr.goes_class,0,1) eq 'M' or strmid(thisflrstr.goes_class,0,1) eq 'X'))
				endif
			endif
		endif

;stop

;skipybin:
		
	endfor

;if i ge 214 then stop

;tvscl,fluximage

;skipxbin:

endfor

if keyword_set(longitude) then filelonlat='_lon' else filelonlat='_lat'
if n_elements(filemod) lt 1 then filemod=''
outfile=datapath+'butterfly_diags'+filelonlat+filemod+'.sav'

save,xarr,yarr,fluximage,areaimage,fluxbalimage,fluxbalunsimage,fluxunipolimage,fluxbalunipolimage,rvalimage,wlsgimage, $
	rvalstarimage, wlsgstarimage, bfluxemergeimage, fluxbipolimage, fluxbalbipolimage, fluxgt22image, fluxgt23image, $
	flrindeximage, flrnumimage, flrmximage, fluxsignedimage, fluxposimage, fluxnegimage, fluxfracimage, fluxfracavgimage, $
	fluxunipolpos, fluxunipolneg, fluxbipolpos, fluxbipolneg, fluxuniemerge, fluxuniemergeavg, file=outfile

return, fluximage


end