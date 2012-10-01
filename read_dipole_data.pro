pro read_dipole_data, out_struct, sav=sav, res1tore=res1tore

datapath='~/science/data/solar_dipole/'
datafile=datapath+'wso_poles_cycle.dat'

if keyword_set(res1tore) then begin
	restore,datapath+'wso_poles.sav'
	return
endif

out_struct={tim:0., bnth:0., bsth:0., bavg:0., filt_nth:0., filt_sth:0., filt_avg:0.}

;readcol,datafile, xx, skip=2, form='A',delim='§‹™'
;readcol,datafile, x1,x2,x3,x4,x5,x6,x7,x8,x9, skip=2, form='A'

;read_ascii, rawstr, file = datafile, format = ['(A22)','(A7)','(A6)','(A8)','(A14)','(A7)','(A7)','(A9)'], skiplines = 2, varnames = varnames, nvar = 8
read_ascii, rawstr, file = datafile, format = '(A80)', skiplines = 2, varnames = varnames, nvar = 1
xx=rawstr.x1

;stop

;1976:05:31_21h:07m:13s    89N -126S  108Avg   20nhz filt:  111Nf  -92Sf  101Avgf
x1=strmid(xx,0,22)
x2=strmid(xx,22,6)
x3=strmid(xx,29,5)
x4=strmid(xx,35,5)
x5=strmid(xx,43,14) ;dummy
x6=strmid(xx,57,5)
x7=strmid(xx,64,5)
x8=strmid(xx,71,5)

yyyy=strmid(x1,0,4)
mo=strmid(x1,5,2)
dd=strmid(x1,8,2)
hh=strmid(x1,11,2)
mm=strmid(x1,15,2)
ss=strmid(x1,19,2)

time=anytim(yyyy+'-'+mo+'-'+dd+'T'+hh+':'+mm+':'+ss+'.000')

;stop

bnth=float(strmid(x2,0,6))
bsth=float(strmid(x3,0,5))
bavg=float(strmid(x4,0,5))

filtnth=float(strmid(x6,0,5))
filtsth=float(strmid(x7,0,5))
filtavg=float(strmid(x7,0,5))

nline=n_elements(x1)

out_struct=replicate(out_struct,nline)
out_struct.tim=time
out_struct.bnth=bnth
out_struct.bsth=bsth
out_struct.bavg=bavg
out_struct.filt_nth=filtnth
out_struct.filt_sth=filtsth
out_struct.filt_avg=filtavg

if keyword_set(sav) then save,out_struct,file=datapath+'wso_poles.sav'

end