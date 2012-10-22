;Inputs: 	xrange = min max time in ???
;			yrange = min max degrees
;			xbin = time in ???
;			ybin = degrees
function smart_bin2d_simple, xarr, yarr, xbin=xbin, ybin=ybin, xrange=xrange, yrange=yrange, years=years, $
	struct=instruct, timarr=timarr, latitude=latitude, longitude=longitude, longarr=longarr, $
	path=path, outfile=outfile,filemod=filemod

if n_elements(path) lt 1 then datapath='~/science/data/cycle23_sav/smart_paper_2/' else datapath=path

struct=instruct

;find number of bins for x and y axes
nx=round(abs(xrange[0]-xrange[1])/float(xbin))
ny=round(abs(yrange[0]-yrange[1])/float(ybin))

;window,xs=nx,ys=ny

fluximage=fltarr(nx,ny)
fluxgt22image=fltarr(nx,ny)
fluxgt5d22image=fltarr(nx,ny)
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
;flrindeximage=fltarr(nx,ny)
;flrnumimage=fltarr(nx,ny)
;flrmximage=fltarr(nx,ny)

if keyword_set(years) then secpyr=3600.*24.*365. else secpyr=1.

if n_elements(timarr) eq 0 then tim=anytim(struct.time)/secpyr else tim=timarr ;-min(anytim(struct.time)/3600./24./365.)
if keyword_set(latitude) then lat=struct.hglat else begin
	if n_elements(longarr) gt 0 then lat=longarr else lat=struct.carlon
endelse

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
				
				;need to calculate each parameter for each unique time and then average over bin. damn it this is not easy.
				utims=thisstruct[uniq(thisstruct.time)].time
				ntims=n_elements(utims)
				thisflux=fltarr(ntims)
				thisarea=fltarr(ntims)
				thisfluxsigned=fltarr(ntims)
				thisfluxpos=fltarr(ntims)
				thisfluxneg=fltarr(ntims)
				thisfluxfrac=fltarr(ntims)
				thisfluxfracavg=fltarr(ntims)
				thisflux22=fltarr(ntims)
				thisflux5d22=fltarr(ntims)
				thisflux23=fltarr(ntims)
				
				;calculate property for each disk image
				for k=0,ntims-1 do begin
					thisthisstr=thisstruct[where(thisstruct.time eq utims[k])]
					
					;Do flux threshold images
					wgt22=where(thisthisstr.bflux ge 1d6)
					wgt5d22=where(thisthisstr.bflux ge 5d6)
					wgt23=where(thisthisstr.bflux ge 1d7)
					if wgt22[0] ne -1 then thisflux22[k]=total(thisthisstr[wgt22].bflux)
					if wgt5d22[0] ne -1 then thisflux5d22[k]=total(thisthisstr[wgt5d22].bflux)
					if wgt23[0] ne -1 then thisflux23[k]=total(thisthisstr[wgt23].bflux)
					
					thisflux[k]=total(thisthisstr.bflux)
					thisarea[k]=total(thisthisstr.area)
					thisfluxsigned[k]=total(thisthisstr.bfluxpos)-abs(total(thisthisstr.bfluxneg))
					thisfluxpos[k]=total(thisthisstr.bfluxpos)
					thisfluxneg[k]=abs(total(thisthisstr.bfluxneg))
					thisfluxfrac[k]=(total(thisthisstr.bfluxpos) < abs(total(thisthisstr.bfluxneg)))/(total(thisthisstr.bfluxpos) > abs(total(thisthisstr.bfluxneg)))
					thisfluxfracavg[k]=mean((thisthisstr.bfluxpos < abs(thisthisstr.bfluxneg))/(thisthisstr.bfluxpos > abs(thisthisstr.bfluxneg)))
				endfor
				
				;average over times
				fluximage[i,j]=mean(thisflux)
				fluxgt22image[i,j]=mean(thisflux22)
				fluxgt5d22image[i,j]=mean(thisflux5d22)
				fluxgt23image[i,j]=mean(thisflux23)
				areaimage[i,j]=mean(thisarea)
				fluxsignedimage[i,j]=mean(thisfluxsigned)
				fluxposimage[i,j]=mean(thisfluxpos)
				fluxnegimage[i,j]=mean(thisfluxneg)
				fluxfracimage[i,j]=mean(thisfluxfrac)
				fluxfracavgimage[i,j]=mean(thisfluxfracavg)
				
	goto,skip_not_properly_binned
				
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
		
	skip_not_properly_binned:
		
			endif
		endif
;skipybin:
		
	endfor

;tvscl,fluximage

;skipxbin:

endfor

if keyword_set(longitude) then filelonlat='_lon' else filelonlat='_lat'
if n_elements(filemod) lt 1 then filemod=''
outfile=datapath+'butterfly_diags'+filelonlat+filemod+'.sav'

save,xarr,yarr,fluximage,areaimage, $;fluxbalimage,fluxbalunsimage,fluxunipolimage,fluxbalunipolimage,rvalimage,wlsgimage, $
	fluxsignedimage, fluxposimage, fluxnegimage, fluxfracimage, fluxfracavgimage, fluxgt22image, fluxgt5d22image, fluxgt23image, $
	file=outfile
	;rvalstarimage, wlsgstarimage, bfluxemergeimage, fluxbipolimage, fluxbalbipolimage, fluxgt22image, fluxgt23image, $
	;flrindeximage, flrnumimage, flrmximage, $
	;fluxunipolpos, fluxunipolneg, fluxbipolpos, fluxbipolneg, fluxuniemerge, fluxuniemergeavg, file=outfile

return, fluximage


end