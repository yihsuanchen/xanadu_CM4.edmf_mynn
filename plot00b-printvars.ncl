load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/ncrc/home2/Yi-hsuan.Chen//script/ncl//ignorantNCL.ncl"
;load "./ignorantNCL.ncl"
load "/ncrc/home2/Yi-hsuan.Chen//script/ncl//ignorantNCL-colormaps.ncl"
;load "zz-read_parameters.ncl"

begin

;------------------
; general settings
;------------------

  wks_name = "test111"
  wks_type = "png"  ; or x11,pdf,ncgm

  datapath = "/lustre/f2/scratch/Yi-hsuan.Chen/work/"

  ;filename = "rfo76b-RF01-mynn_DIFF_M3_EDMFimpCD_NOsm_1x0m1d_1x1a.o268893242/output.stager/lustre/f2/dev/Yi-hsuan.Chen/awg/xanadu/SCM_am4_xanadu_edmf_mynn/rfo76b-RF01-mynn_DIFF_M3_EDMFimpCD_NOsm/ncrc4.intel19-prod-openmp/1x0m1d_1x1a/history/2001071000.nc/20010710.atmos_edmf_mynn.nc"

  filename = "b70a-BOMEX-am4p0_1x0m1d_1x1a.o268896217/output.stager/lustre/f2/dev/Yi-hsuan.Chen/awg/xanadu/SCM_am4_xanadu_edmf_mynn/b70a-BOMEX-am4p0/ncrc4.intel19-prod-openmp/1x0m1d_1x1a/history/1969062200.nc/19690622.atmos_uw.nc"

  infilename = datapath+"/"+filename

  infile = addfile(infilename, "r")


  ;choice = "qadt_edmf_MF_all"
  ;choice = "qldt_edmf_MF_all"
  ;choice = "edmf_a_w_qc"  
  choice = "qcdt_uw"

  choice = "cldliq"

  ;choice = "varname"
  ;varname = "liq_wat"

  if (choice.eq."xlklskdflksld") then
    aa=0

  elseif (choice.eq."cldliq") then
    liq_wat = infile->qldt_uws
    qadt_uws = infile->qadt_uws
    qndt_uws = infile->qndt_uws

    print(choice+", qldt_uws: "+qldt_uws+", qadt_uws: "+qadt_uws+", qndt_uws:"+qndt_uws)

  elseif (choice.eq."qcdt_uw") then
    qldt_uws = infile->qldt_uws
    qadt_uws = infile->qadt_uws
    qndt_uws = infile->qndt_uws

    print(choice+", qldt_uws: "+qldt_uws+", qadt_uws: "+qadt_uws+", qndt_uws:"+qndt_uws)

  elseif (choice.eq."qadt_edmf_MF_all") then
    qadt_edmf_MF = infile->qadt_edmf_MF
    qadt_edmf_MF_adv = infile->qadt_edmf_MF_adv
    qadt_edmf_MF_ent = infile->qadt_edmf_MF_ent
    qadt_edmf_MF_eddy = infile->qadt_edmf_MF_eddy
    qadt_edmf_MF_advPeddy = plus(qadt_edmf_MF_adv, qadt_edmf_MF_eddy)
    qadt_edmf_MF_all  = plus(qadt_edmf_MF_advPeddy, qadt_edmf_MF_ent)
    print(choice+", MF:"+qadt_edmf_MF+", all: "+qadt_edmf_MF_all+" , adv:"+qadt_edmf_MF_adv+", eddy:"+qadt_edmf_MF_eddy+", ent:"+qadt_edmf_MF_ent)

  elseif (choice.eq."qldt_edmf_MF_all") then
    qldt_edmf_MF = infile->qldt_edmf_MF
    qldt_edmf_MF_adv = infile->qldt_edmf_MF_adv
    qldt_edmf_MF_ent = infile->qldt_edmf_MF_ent
    qldt_edmf_MF_eddy = infile->qldt_edmf_MF_eddy
    qldt_edmf_MF_advPeddy = plus(qldt_edmf_MF_adv, qldt_edmf_MF_eddy)
    qldt_edmf_MF_all  = plus(qldt_edmf_MF_advPeddy, qldt_edmf_MF_ent)
    print(choice+", MF:"+qldt_edmf_MF+", all: "+qldt_edmf_MF_all+" , adv:"+qldt_edmf_MF_adv+", eddy:"+qldt_edmf_MF_eddy+", ent:"+qldt_edmf_MF_ent)

  elseif (choice.eq."qadt_edmf_ALL") then
    qadt_edmf = infile->qadt_edmf
    qadt_edmf_ED = infile->qadt_edmf_ED
    qadt_edmf_MF = infile->qadt_edmf_MF
    print(choice+", EDMF:"+qadt_edmf+" , ED:"+qadt_edmf_ED+", MF:"+qadt_edmf_MF)

  elseif (choice.eq."edmf_a_w_qc") then
    edmf_a  = infile->edmf_a
    edmf_w  = infile->edmf_w
    edmf_qc = infile->edmf_qc
    print(choice+", edmf_a:"+edmf_a+", edmf_w: "+edmf_w+", edmf_qc: "+edmf_qc)

  elseif (choice.eq."varname") then
   var1 = infile->$varname$
   print(var1)

  end if

end
