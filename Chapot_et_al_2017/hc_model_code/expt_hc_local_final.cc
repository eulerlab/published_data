/* Experiment cbp_flash */
/*  for nc script retsim.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

int conearr;
int scone_id;
int elnode;

double temp_freq;
double ntrials;
double dstim;
double sdia;
double stimtime;
double inten;
double inten2;
double minten;
double setploti;
double tau;
double tstep;
int nrepeats;
int stimtype;
double predur;
double pulsedur;
double pulsegap;
double prestimdur;
double stimdur;
double poststimdur;
double start_off;
double s_inten;
double m_inten;
double kdist, kprox, ksoma, cadist, caprox, casoma, capdt, cappx, capsm, capkdt, capkpx, capksm;
double dvshb, dvsc, dcrm, drm;

double istart;
double istop;
double istep;
double ipre;
double elec_rs;
double elec_cap;
int stimtip;

int rectype;


/*------------------------------------------------------*/

void defparams(void) 
{
  setptr("conearr",   &conearr);
  setptr("scone_id",  &scone_id);
  
  setptr("temp_freq", &temp_freq);
  setptr("tstep",     &tstep);
  setptr("tau",       &tau);
  setptr("ntrials",   &ntrials);
  setptr("dstim",     &dstim);
  setptr("sdia",      &sdia);
  setptr("stimtime",  &stimtime);
  setptr("minten",    &minten);
  setptr("setploti",  &setploti);
  setptr("nrepeats",  &nrepeats);
  setptr("stimtype",  &stimtype);
  setptr("predur",    &predur);
  setptr("pulsedur",  &pulsedur);
  setptr("pulsegap",  &pulsegap);
  setptr("start_off", &start_off);
  setptr("s_inten",   &s_inten);
  setptr("m_inten",   &m_inten);
  setptr("kdist",     &kdist);
  setptr("kprox",     &kprox);
  setptr("ksoma",     &ksoma);
  setptr("cadist",    &cadist);
  setptr("caprox",    &caprox);
  setptr("casoma",    &casoma);
  setptr("capdt",     &capdt);
  setptr("cappx",     &cappx);
  setptr("capsm",     &capsm);
  setptr("capkdt",    &capkdt);
  setptr("capkpx",    &capkpx);
  setptr("capksm",    &capksm);
  setptr("dvshb",     &dvshb);
  setptr("dvsc",      &dvsc);
  setptr("dcrm",      &dcrm);
  setptr("drm",       &drm);

  setptr("predur",	&predur);
  setptr("prestimdur",	&prestimdur);
  setptr("stimdur",	&stimdur);
  setptr("poststimdur",	&poststimdur);
  setptr("istart",	&istart);
  setptr("istop",	&istop);
  setptr("istep",	&istep);
  setptr("ipre",	&ipre);
  setptr("stimtip",	&stimtip);

  setptr("rectype",&rectype);
  
  nvalfile = "nval_hb_final.n";
  hb_file="morph_hc_dendrite";
  hb_densfile="dens_hb_final.n";
  hb_densfile2="dens_hb_final.n";
  cone_densfile="dens_cone.n";
  
}

/*------------------------------------------------------*/

void setparams(void)
{
  make_rods = 0;
  if (stimtype != 2) make_cones= 1;        /* make cones, dbp, gc */
  else make_cones=0;
  make_ha   = 0;
  make_hb   = 1;
  make_hbat = 0;
  make_dbp1 = 0;
  make_dbp2 = 0;
  make_rbp  = 0;
  make_sbac = 0;
  make_gca  = 0;
  make_gcb  = 0;
  make_dsgc = 0;

  make_cone_hb = 0;
  make_hb_hb   = 0;
  n_cones = 1; //one dummy cone needed to prevent make_cones=0;

  if (notinit(bg_inten)) bg_inten = 2.5e3;
  if (notinit(kdist))       kdist = 0.0e-3;
  if (notinit(kprox))       kprox = 0.0e-3;
  if (notinit(ksoma))       ksoma = 0.0e-3;
  if (notinit(cadist))     cadist = 0.0e-3;
  if (notinit(caprox))     caprox = 0.0e-3;
  if (notinit(casoma))     casoma = 0.0e-3;
  if (notinit(capdt))       capdt = 0.0e-7;
  if (notinit(cappx))       cappx = 0.0e-7;
  if (notinit(capsm))       capsm = 0.0e-7;
  if (notinit(capkdt))     capkdt = 0.0e-7;
  if (notinit(capkpx))     capkpx = 0.0e-7;
  if (notinit(capksm))     capksm = 0.0e-7;

  //Set Cone and Horizontal Cell Resting Potentials
  
  vk = -0.082;
  if (notinit(dvsc))  dvsc  = -0.04;		// for dens_cone.n
  if (notinit(dcrm))   dcrm = 1.2e3;
  if (notinit(dvshb)) dvshb = -0.04;		// for dens_hb_final.n
  if (notinit(drm))     drm = 2.5e3;		 
  if (notinit(dri))     dri = 200;

  setn(hb,NCOLOR,RCOLOR);

}

/*------------------------------------------------------*/

void setdens(void)
/* set density and conductance parameters */

{

}

/*------------------------------------------------------*/

void addcells(void)
{ 
  if (stimtype!=2) {
  checkcellout(xcone); //remove dummy cone
  if (notinit (n_cones)) n_cones = 1;
  if (notinit(conearr)) conearr = 1;
  if (notinit(scone_id)) scone_id = -1;
  node *npnt;
  int position_nodes[]={94,682,276,353,401,457,511,651,764,788};
  int contact_cones[]={1,2,3,4,5,6,6,7,8,2,2,7,9,10,1,1};
  int contact_nodes[]={94,120,276,353,401,457,474,511,651,679,682,696,764,788,817,819};

  #define CONEARR 30
  if (!notinit(conearr)) {
    conexarr  = (double *)emalloc(CONEARR*sizeof(double));
    coneyarr  = (double *)emalloc(CONEARR*sizeof(double));
    conetharr = (double *)emalloc(CONEARR*sizeof(double));
    conenarr  = (int *) emalloc(CONEARR*sizeof(int));
    if (conearr==0) {
        n_cones = 0;
    }
    if (conearr==1) {
      n_cones=10;
      for(int i=0;i<n_cones;i++) {
	if(i!=scone_id-1) cone_pigm=11;
	else cone_pigm=12;
	npnt=ndn(hb,1,position_nodes[i]);
        conexarr[i]=npnt->xloc; coneyarr[i]=npnt->yloc; conetharr[i]=0; conenarr[i]=i+1;
	ncones+=setupcells(xcone,1,conexarr+i,coneyarr+i,conetharr+i,conenarr+i);
      }
    }
  }
 
  fprintf(stderr,"ncones: %i, n_cones: %i\n",ncones,n_cones);  

  for(int i=0;i<16;i++) {
    connect_synapse(xcone,contact_cones[i],1,hb,1,contact_nodes[i]);
  }
  }
}

/*------------------------------------------------------*/

void runexpt(void)

{
    int ct, cn, n, r, plnum, electrode_node;
    int colr;
    int synin1, synin2;
    double t, fmax,fmin, tim;
    double m, rmin, rmax, plsize;
    double dtrial,exptdur;
    double Vmin, Vmax;
    double maxca;
    double Imin, Imax;
    double ipulse;
    int s_wavel=430;
    int m_wavel=527;
#define GAUSRNG 1

  
  if (notinit(temp_freq)) temp_freq = 2;

  if (temp_freq == 0) {
    fprintf (stderr,"## retsim1: temp_freq = 0, STOPPING\n");
    temp_freq = 1;
  };

  ploti = 1e-3;
  electrode_node = 5000;

  if (notinit(stimtype))   stimtype = 1;
  if (notinit(tstep))         tstep = 0.001;     /* stimulus sampling period */
  if (notinit(tau))             tau = 0.01;      /* stimulus lowpass filter period */
  if (notinit(exptdur))     exptdur = 5;         /* expt duration */
  if (notinit(sdia))           sdia = 500;       /* spot diameter */
  if (notinit(stimtime))   stimtime = 0.0;       /* stimulus time */
  if (notinit(minten))       minten = bg_inten;  /* background intensity (for make cone)*/
  if (!notinit(setploti))     ploti = setploti;  /* plot time increment */
  if (stimtype == 2) if (notinit(nrepeats))   nrepeats = 1;         /* number of stimulus repeats */
  else if (notinit(nrepeats))   nrepeats = 3;
  if (notinit(pulsedur))   pulsedur = 1;
  if (notinit(dtrial))       dtrial = 15;
  if (notinit(pulsegap))   pulsegap = 4;
  if (notinit(start_off)) start_off = 0.5;
  if (notinit(s_inten))     s_inten = 13;
  if (notinit(m_inten))     m_inten = 13;
  if (notinit(rectype))   rectype = 0;

  if (notinit(ipre))          ipre  = 0e-12;
  if (notinit(istart))       istart =   20e-13;
  if (notinit(istop))         istop =   100e-13;
  if (notinit(istep))         istep =   20e-13;
  if (notinit(stimtip))  stimtip  = 2;


  if (dtrial<pulsedur) dtrial=pulsedur;

  endexp  = exptdur*nrepeats;


  if (stimtype!=2) {
    plot_v_nod(ct=xcone,cn=5,n=1,Vmin=-.05,Vmax=-.02,colr=blue,"",-1,0.5);
    if (scone_id>0) plot_v_nod(ct=xcone,cn=scone_id,n=1,Vmin=-.05,Vmax=-.02,colr=blue,"",-1,0.5);
    plot_ca_nod(ct=xcone,cn=5,n=1,maxca=20e-4,colr=brown,"",5,0.5);
    if (scone_id>0) plot_ca_nod(ct=xcone,cn=scone_id,n=1,maxca=20e-4,colr=brown,"",5,0.5);
    }
  
  int pn0[]={94,120,276,353,401,457,474,511,651,679,682,696,764,788,817,819,0};
  int pn2[]={682,101,53,45,37,29,18,10,4,0};
  int pn3[]={276,267,212,203,193,185,175,166,157,29,18,10,4,0};
  int pn20[]={0,94,682,276,353,401,457,511,651,764,788};
  int pn100[824];
  for(int i=0;i<824;i++) pn100[i]=i;
  int *plot_nodes;
  int pn_size=0;
  plot_nodes=0;
	  
  switch(rectype) {
    case 0:
      plot_nodes=pn0;
      pn_size=sizeof(pn0) / sizeof(pn0[0]);
      break;
    case 2:
      plot_nodes=pn2;
      pn_size=sizeof(pn2) / sizeof(pn2[0]);
      break;
    case 3:
      plot_nodes=pn3;
      pn_size=sizeof(pn3) / sizeof(pn3[0]);
      break;
    case 30:
      plot_nodes=pn20;
      pn_size=sizeof(pn20) / sizeof(pn20[0]);
      break;
    case 20:
      plot_nodes=&pn20[stimtip];
      pn_size=1;
      break;
    case 100:
      plot_nodes=pn100;
      pn_size=100;
      break;
    case 101:
      plot_nodes=&pn100[100];
      pn_size=100;
      break;
    case 102:
      plot_nodes=&pn100[200];
      pn_size=100;
      break;
    case 103:
      plot_nodes=&pn100[300];
      pn_size=100;
      break;
    case 104:
      plot_nodes=&pn100[400];
      pn_size=100;
      break;
    case 105:
      plot_nodes=&pn100[500];
      pn_size=100;
      break;
    case 106:
      plot_nodes=&pn100[600];
      pn_size=100;
      break;
    case 107:
      plot_nodes=&pn100[700];
      pn_size=100;
      break;
    case 108:
      plot_nodes=&pn100[800];
      pn_size=24;
      break;
  }

  for (int i=0;i<pn_size;i++) plot_v_nod(ct=hb,cn=1,n=plot_nodes[i],Vmin=-.05,Vmax=-.02,colr=red,"",-1,0.5);
  for (int i=0;i<pn_size;i++) plot_ca_nod(ct=hb,cn=1,n=plot_nodes[i],1,maxca=20e-4,colr=green,"",5,0.5);

  if (stimtype == 1) {
    stim_backgr(minten);
    dtrial=3*(pulsedur+pulsegap);
    exptdur=nrepeats*dtrial;
    for (t=n=0; n<nrepeats; n++,t+=dtrial){
      double start, dur;
      stim_spot(sdia, 0, 0, m_inten*minten, start=t+stimtime+start_off,dur=pulsedur,m_wavel,0); 
      stim_spot(sdia, 0, 0, m_inten*minten, start=t+stimtime+start_off+pulsedur+pulsegap,dur=pulsedur,s_wavel,0.0);
      stim_spot(sdia, 0, 0, m_inten*minten, start=t+stimtime+start_off+2*(pulsedur+pulsegap),dur=pulsedur,m_wavel,0); 
      stim_spot(sdia, 0, 0, m_inten*minten, start=t+stimtime+start_off+2*(pulsedur+pulsegap),dur=pulsedur,s_wavel,0.0);	
      step(dtrial);
    }
  }


  else if (stimtype == 2) { //current clamp
    ct=hb;
    cn=1;
    if (notinit(prestimdur))   prestimdur = 0.05;
    if (notinit(stimdur))      stimdur    = 0.2;
    if (notinit(poststimdur)) poststimdur = 0.5;

    if (notinit(elec_rs))  elec_rs   = 20e6;
    if (notinit(elec_cap)) elec_cap  = 1e-14;
    if (notinit(elnode))     elnode  = soma;

    int tip_nodes[]={0,94,682,276,353,401,457,511,651,764,788};
    elnode = tip_nodes[stimtip];

    if (notinit(predur)) predur = 0.15;

    simtime = 0 - predur;
    step (predur); 
    int i;

    for (i=0,ipulse=istart; ipulse <= istop; i++,ipulse += istep) {
      simtime = 0;
      cclamp(ndn(ct,cn,elnode), ipre,   simtime,  prestimdur);
      step (prestimdur);

      cclamp(ndn(ct,cn,elnode), ipulse, simtime,  stimdur);
      step (stimdur);

      cclamp(ndn(ct,cn,elnode), ipre, simtime,  poststimdur);
      step (poststimdur); 
    }
  }
}
