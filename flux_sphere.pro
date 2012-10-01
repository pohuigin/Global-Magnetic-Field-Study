;measure the flux on a portion of a sphere with some heliographic bounds
;BAVG=average befield over region (in Gauss)
;	1 element or qd array
;LONRANGE=min and max heliographic longitudes of region (in degrees)
;	2 elements or [2,*] element array
;LATRANGE=min and max heliographic latitudes of region (in degrees)
;	2 elements or [2,*] element array
;
;Output: Flux in Mx (cm^2 G)
;	1 elements or 1d array

function flux_sphere,bavg,lonrange=lonrange,latrange=latrange

rsun_Mm=6.955d10 ;in cm

if n_elements(size(lonrange,/dim)) eq 1 or n_elements(size(latrange,/dim)) eq 1 then begin
	phi_ba=[lonrange[1],lonrange[0]]*!dtor
	theta_ba=[latrange[1],latrange[0]]*!dtor
	flux=-(phi_ba[0]-phi_ba[1])*(cos(theta_ba[0])-cos(theta_ba[1]))*bavg*rsun_Mm^(2.)
endif else begin
	phi_ba=[lonrange[1,*],lonrange[0,*]]*!dtor
	theta_ba=[latrange[1,*],latrange[0,*]]*!dtor
	flux=-(phi_ba[0,*]-phi_ba[1,*])*(cos(theta_ba[0,*])-cos(theta_ba[1,*]))*bavg*rsun_Mm^(2.)
endelse



return,flux

end