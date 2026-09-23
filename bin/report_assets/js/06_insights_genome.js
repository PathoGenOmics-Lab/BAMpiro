function insGenome(){
  var S=((R.samples)||[]).filter(onGenome); if(!S.length) return '';   // the reference in view, where there are several
  var thr=R.thresholds||{};
  var G=genome(), gl=G.len||0;
  // ---- cohort callability: prefer callable_pct, else 100 - missing_pct ----
  var cvals=[];
  for(var i=0;i<S.length;i++){var m=S[i].m||{},c=null;
    if(m.callable_pct!=null)c=m.callable_pct;
    else if(m.missing_pct!=null)c=100-m.missing_pct;
    if(c!=null)cvals.push(c);}
  var medCall=cvals.length?_median(cvals):null;
  // ---- low-callability windows: cohort mean missing per bin vs missing_max cut-off (masked bins excluded, mirroring the plot) ----
  var missRows=S.filter(function(s){return s.miss&&s.miss.length;});
  var nb=missRows.length?missRows[0].miss.length:0;
  var mb=(G.mask_bins&&G.mask_bins.length===nb)?G.mask_bins:null;
  var missMaxPct=(thr.missing_max!=null)?thr.missing_max:10;
  var lowN=0,evalBins=0;
  for(var b=0;b<nb;b++){
    if(mb&&mb[b]>=0.5)continue;
    var sum=0,cnt=0;
    for(var r=0;r<missRows.length;r++){var vv=missRows[r].miss[b];if(vv!=null){sum+=vv;cnt++;}}
    if(!cnt)continue; evalBins++;
    if((sum/cnt)*100>missMaxPct)lowN++;
  }
  // ---- densest SNP window across the cohort -> gene / region ----
  var snpRows=S.filter(function(s){return s.trk&&s.trk.snp&&s.trk.snp.length;});
  var snb=snpRows.length?snpRows[0].trk.snp.length:0;
  var peakBin=-1,peakVal=-1;
  for(var b2=0;b2<snb;b2++){var s2=0,c2=0;
    for(var r2=0;r2<snpRows.length;r2++){var v2=snpRows[r2].trk.snp[b2];if(v2!=null){s2+=v2;c2++;}}
    if(!c2)continue; var avg=s2/c2; if(avg>peakVal){peakVal=avg;peakBin=b2;}}
  var peakGene=null,pMb=0,p0Mb=0,p1Mb=0;
  if(peakBin>=0&&snb&&gl){
    var p0=peakBin/snb*gl,p1=(peakBin+1)/snb*gl; pMb=(p0+p1)/2/1e6; p0Mb=p0/1e6; p1Mb=p1/1e6;
    var genes=G.genes||[];
    for(var g=0;g<genes.length;g++){var ge=genes[g];
      if(ge.start!=null&&ge.end!=null&&ge.start<=p1&&ge.end>=p0){peakGene=ge.name;break;}}
  }
  if(medCall==null&&!evalBins&&peakBin<0)return '';
  var peakR=Math.round(peakVal*10)/10;
  // ---- narrative ----
  var ctone=(medCall==null)?'':(medCall<80?' class="tone-bad"':medCall<90?' class="tone-warn"':'');
  var narr='';
  if(medCall!=null){
    narr+='Cohort median callability is <b'+ctone+'>'+medCall.toFixed(1)+'%</b>'+(gl?(' of the '+(gl/1e6).toFixed(2)+' Mb reference'+(genomeRefs().length>1?(' '+esc(genomeRef())):'')):'')+'.';
  }
  if(evalBins){
    narr+=' '+(lowN>0
      ?('<b class="tone-warn">'+lowN+' of '+evalBins+'</b> windows sit below the &gt;'+missMaxPct+'% missing cut-off')
      :('all '+evalBins+' scored windows clear the &gt;'+missMaxPct+'% missing cut-off'))+'.';
  }
  if(peakBin>=0){
    narr+=' Cohort SNP density peaks '+(peakGene?('in <b>'+esc(peakGene)+'</b>'):('in the '+p0Mb.toFixed(2)+'–'+p1Mb.toFixed(2)+' Mb window'))+' (~'+pMb.toFixed(2)+' Mb, '+peakR+' SNPs/bin) — a coarse '+(gl&&snb?(Math.round(gl/snb)+' bp'):'')+' cohort mean, not a per-site call.';
  }
  // ---- chips ----
  var chips=[];
  if(medCall!=null)chips.push({t:medCall.toFixed(1)+'% callable',cls:medCall<80?'bad':medCall<90?'warn':'good',title:'cohort median of per-sample callable fraction'});
  if(evalBins)chips.push({t:lowN+' low-call window'+(lowN===1?'':'s'),cls:lowN>0?'warn':'good',title:'cohort mean N-calls > '+missMaxPct+'% across '+evalBins+' scored bins (masked bins excluded)'});
  if(peakBin>=0)chips.push({t:'SNP peak '+(peakGene?esc(peakGene):(pMb.toFixed(2)+' Mb')),cls:'',title:peakR+' SNPs/bin, cohort mean over ~'+(gl&&snb?Math.round(gl/snb):0)+' bp bins'});
  return insBox('Callability',narr,chips);
}
function insFunction(){
  var rows=R.samples.filter(function(s){return s.m&&s.m.ann_high!=null;});
  if(!rows.length)return "";
  var totHigh=0;rows.forEach(function(s){totHigh+=(s.m.ann_high||0);});
  var lofVals=[];R.samples.forEach(function(s){if(s.m&&s.m.lof_pct!=null)lofVals.push(s.m.lof_pct);});
  var med=lofVals.length?_median(lofVals):null;
  var cut=(med!=null)?Math.max(1,2*med):null;
  var hot=[];
  if(cut!=null){
    R.samples.forEach(function(s){var v=(s.m&&s.m.lof_pct!=null)?s.m.lof_pct:null;
      if(v!=null&&v>0&&v>=cut)hot.push({s:s.s,v:v});});
    hot.sort(function(a,b){return b.v-a.v;});
  }
  var narr='Cohort carries <b>'+Math.round(totHigh).toLocaleString('en-US')+'</b> HIGH-impact variant'+(Math.round(totHigh)==1?'':'s')+' across <b>'+rows.length+'</b> annotated sample'+(rows.length==1?'':'s')+' (snpEff)';
  if(med!=null)narr+=', median LOF fraction <b>'+med.toFixed(1)+'%</b>';
  narr+='. ';
  var chips=[];
  if(hot.length){
    narr+='<b class="tone-warn">'+hot.length+'</b> sample'+(hot.length==1?'':'s')+' show elevated loss-of-function (&ge;'+cut.toFixed(1)+'%, &gt;2&times; median) &mdash; watch for pseudogenisation or indel/annotation artifact.';
    hot.slice(0,6).forEach(function(h){
      chips.push({t:h.s+' '+h.v.toFixed(1)+'%',cls:(h.v>=Math.max(2,4*(med||0))?'bad':'warn'),title:'lof_pct '+h.v.toFixed(2)+'% (cohort median '+((med||0)).toFixed(2)+'%)'});
    });
    if(hot.length>6)chips.push({t:'+'+(hot.length-6)+' more',cls:'warn',title:'additional elevated-LOF samples'});
  } else {
    narr+='LOF fraction is uniform across samples &mdash; no pseudogenisation outlier stands out.';
  }
  return insBox("Impact",narr,chips.length?chips:null);
}
function insGeneBurden(){
  var gb=R.gene_burden; if(!(gb&&gb.length))return '';
  var scored=[],i,g,burd;
  for(i=0;i<gb.length;i++){g=gb[i];
    burd=(g.total_impactful!=null?g.total_impactful:(g.high||0)+(g.moderate||0));
    if(burd>0)scored.push({gene:g.gene,burd:burd,high:(g.high||0),mod:(g.moderate||0),ns:(g.n_samples==null?null:g.n_samples),eff:g.dominant_effect||'NA'});
  }
  if(!scored.length)return '';
  scored.sort(function(a,b){return (b.burd-a.burd)||(b.high-a.high);});
  var top=scored.slice(0,3);
  var ncoh=(R.samples&&R.samples.length)||0;
  var totalHigh=0; for(i=0;i<scored.length;i++)totalHigh+=scored[i].high;
  var lead=top[0];
  var names=top.map(function(t){return '<b>'+esc(t.gene)+'</b>';});
  var namestr=(names.length===1?names[0]:(names.slice(0,-1).join(', ')+' and '+names[names.length-1]));
  var narr='Impactful (HIGH+MODERATE) variants concentrate in '+namestr+'. '+
    '<b>'+esc(lead.gene)+'</b> leads with '+lead.burd+' impactful variant'+(lead.burd===1?'':'s')+
    (lead.ns!=null?(' across '+lead.ns+(ncoh?(' of '+ncoh):'')+' sample'+(lead.ns===1?'':'s')):'')+
    (lead.eff&&lead.eff!=='NA'?(' (dominant effect '+esc(lead.eff)+')'):'')+'.';
  if(totalHigh>0){
    var hg=[]; for(i=0;i<scored.length;i++)if(scored[i].high>0)hg.push(scored[i].gene);
    narr+=' <b class="tone-bad">'+totalHigh+' HIGH-impact</b> (predicted loss-of-function) call'+(totalHigh===1?'':'s')+' present';
    if(hg.length)narr+=' in '+hg.map(function(x){return '<b>'+esc(x)+'</b>';}).join(', ');
    narr+='.';
  } else {
    narr+=' No HIGH-impact calls — burden is dominated by missense / in-frame changes.';
  }
  narr+=' This counts variant sites per gene, not per-sample QC flags.';
  var chips=top.map(function(t){
    return {t:t.gene+' · '+t.burd+' var'+(t.ns!=null?(', '+t.ns+' smp'):''),
            cls:(t.high>0?'bad':'warn'),
            title:esc(t.gene)+': '+t.high+' HIGH, '+t.mod+' MODERATE'+(t.ns!=null?(', '+t.ns+' samples'):'')+', dominant '+esc(t.eff)};
  });
  return insBox('Gene burden', narr, chips);
}
function insHotspots(){
  var S=(R.samples||[]).filter(onGenome); if(!S.length) return '';
  var agg=null,ns=0;
  for(var i=0;i<S.length;i++){
    var t=S[i].trk&&S[i].trk.snp; if(!t||!t.length) continue;
    if(!agg){agg=[];for(var k=0;k<t.length;k++)agg.push(0);}
    var n=Math.min(agg.length,t.length);
    for(var j=0;j<n;j++)agg[j]+=(+t[j]||0);
    ns++;
  }
  if(!agg||!ns) return '';
  var G=genome(), nb=R.nbins||agg.length, gl=G.len||nb, binbp=gl/nb;
  function b0of(p){return Math.max(0,Math.min(nb-1,Math.floor(p/gl*nb)));}
  var genes=G.genes||[]; if(!genes.length) return '';
  var items=[];
  for(var gi=0;gi<genes.length;gi++){
    var g=genes[gi],s=Math.max(1,g.start),e=Math.max(s,g.end),len=(e-s+1)||1,b0=b0of(s),b1=b0of(e),tot=0;
    for(var b=b0;b<=b1;b++){var binS=b*binbp,binE=(b+1)*binbp,ov=Math.max(0,Math.min(e,binE)-Math.max(s,binS));tot+=(agg[b]||0)*(binbp>0?ov/binbp:1);}
    items.push({name:g.name,snp:tot,dens:tot/(len/1000)});
  }
  if(!items.length) return '';
  var dv=items.map(function(x){return x.dens;}),med=_median(dv),ad=[];
  for(var d=0;d<dv.length;d++)ad.push(Math.abs(dv[d]-med));
  var spread=_median(ad)*1.4826;
  items.forEach(function(x){x.z=spread?(x.dens-med)/spread:0;});
  items.sort(function(a,b){return b.dens-a.dens;});
  var top=items[0]; if(!(top.dens>0)) return '';
  var out=items.filter(function(x){return x.z>=3;});
  var tone=(top.z>=3)?'tone-bad':(top.z>=2?'tone-warn':'');
  var zs=(top.z>0?'+':'')+top.z.toFixed(1);
  var narr='Cohort SNP density peaks at <b'+(tone?' class="'+tone+'"':'')+'>'+esc(top.name)+'</b> ('+top.dens.toFixed(1)+' SNPs/kb, robust z '+zs+' vs the gene distribution). '+
    (out.length>1?('<b class="tone-warn">'+out.length+'</b> genes clear z=3. '):'')+
    'A density spike can be positive selection, recombination, or a masking / mapping artifact - confirm before reading it as biology.';
  var chips=[],picks=out.length?out:[top];
  for(var p=0;p<picks.length&&p<4;p++){
    chips.push({t:picks[p].name+' '+picks[p].dens.toFixed(1)+'/kb',cls:picks[p].z>=3?'bad':(picks[p].z>=2?'warn':''),title:'robust z '+(picks[p].z>0?'+':'')+picks[p].z.toFixed(1)+' over '+ns+' samples'});
  }
  return insBox('Hotspots',narr,chips);
}
function insTemporal(){
  var S=(R&&R.samples)||[];
  function yrOf(d){if(d==null)return null;var m=String(d).match(/\d{4}/);if(!m)return null;var y=+m[0];return (y>=1000&&y<=2100)?y:null;}
  var years=[],pts=[],i,s,y,sn;
  for(i=0;i<S.length;i++){s=S[i];y=yrOf(s.date);if(y==null)continue;years.push(y);sn=(s.m&&s.m.snps!=null)?s.m.snps:null;if(sn!=null)pts.push({y:y,v:sn});}
  if(!years.length)return '';
  var ymin=Math.min.apply(null,years),ymax=Math.max.apply(null,years),span=ymax-ymin;
  var nDated=years.length,M=S.length;
  // cross-sectional OLS slope of SNP proxy vs collection year (only if enough spread)
  var uy={},k;for(k=0;k<pts.length;k++){uy[pts[k].y]=1;}
  var distinctY=0;for(k in uy){if(uy.hasOwnProperty(k))distinctY++;}
  var slope=null,r=null;
  if(pts.length>=3&&distinctY>=2){
    var n=pts.length,sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(k=0;k<n;k++){var px=pts[k].y,pv=pts[k].v;sx+=px;sy+=pv;sxx+=px*px;syy+=pv*pv;sxy+=px*pv;}
    var den=n*sxx-sx*sx;
    if(den>0){slope=(n*sxy-sx*sy)/den;var rd=Math.sqrt(den*(n*syy-sy*sy));if(rd>0)r=(n*sxy-sx*sy)/rd;}
  }
  var narr='Cohort spans <b>'+span+' yr</b> ('+ymin+'–'+ymax+'); '+nDated+' of '+M+' samples carry a parseable collection year';
  if(span===0){
    narr+=', all from a <b class="tone-warn">single year</b> — no temporal spread to resolve an accumulation trend.';
  }else if(slope!=null){
    var per=(slope>=0?'+':'')+slope.toFixed(1);
    narr+='. Cross-sectionally the SNP proxy tracks <b>'+per+' SNPs/yr</b>'+(r!=null?' (r='+r.toFixed(2)+', n='+pts.length+')':'')+' — a rough correlation, <b class="tone-warn">not a molecular clock</b>: it mixes lineages and is unrooted, so read it as a sanity check, not a rate.';
  }else{
    narr+='; <b class="tone-warn">too few dated points</b> ('+pts.length+') to fit a SNP-accumulation trend.';
  }
  var chips=[{t:span+' yr span',cls:'',title:ymin+'–'+ymax+' collection years'},{t:nDated+'/'+M+' dated',cls:(nDated<M?'warn':'good'),title:'samples with a parseable collection year'}];
  if(slope!=null)chips.push({t:(slope>=0?'+':'')+slope.toFixed(1)+' SNP/yr',cls:'',title:'cross-sectional OLS slope of SNP proxy vs year'+(r!=null?', r='+r.toFixed(2):'')+' (heuristic, not a clock)'});
  return insBox('Temporal',narr,chips);
}
function insPnps(){
  var pp=R.pnps; if(!(pp&&pp.length))return '';
  var pos=[],pur=[],neutral=0,scored=0,i;
  for(i=0;i<pp.length;i++){
    var g=pp[i]; if(g.pnps==null)continue; scored++;
    if(g.pnps>1)pos.push(g);
    else if(g.pnps<1)pur.push(g);
    else neutral++;
  }
  if(!scored)return '';
  pos.sort(function(a,b){return b.pnps-a.pnps;});
  pur.sort(function(a,b){return a.pnps-b.pnps;});
  var chips=[],narr,j;
  if(pos.length){
    var names=pos.slice(0,3).map(function(x){return esc(x.gene);}).join(', ');
    var top=pos[0];
    narr='<b class="tone-warn">'+pos.length+'</b> of '+scored+' scored gene'+(scored===1?'':'s')+' show pN/pS &gt; 1 ('+names+(pos.length>3?', +'+(pos.length-3)+' more':'')+'), a signal of <b>candidate positive selection</b> - <b>'+esc(top.gene)+'</b> highest at '+top.pnps.toFixed(2)+(pur.length?'; the other '+pur.length+' skew purifying':'')+'. Read as a within-sample proxy, <b>not</b> a phylogenetic dN/dS.';
    for(j=0;j<pos.length&&j<6;j++){
      chips.push({t:pos[j].gene+' '+pos[j].pnps.toFixed(2),cls:'warn',title:'pN='+(pos[j].pn==null?'NA':pos[j].pn.toFixed(3))+' pS='+(pos[j].ps==null?'NA':pos[j].ps.toFixed(3))+(pos[j].n==null?'':' over '+Math.round(pos[j].n)+' pairs')+' - '+(pos[j].effect||'NA')});
    }
  } else {
    var lowest=pur.length?pur[0]:null;
    narr='All '+scored+' scored gene'+(scored===1?'':'s')+' sit at pN/pS &le; 1'+(lowest?', lowest <b>'+esc(lowest.gene)+'</b> at '+lowest.pnps.toFixed(2):'')+' - dominated by <b>purifying selection</b>, no positive-selection candidates. Within-sample proxy, <b>not</b> a phylogenetic dN/dS.';
    for(j=0;j<pur.length&&j<6;j++){
      chips.push({t:pur[j].gene+' '+pur[j].pnps.toFixed(2),cls:'good',title:'purifying'+(pur[j].n==null?'':' over '+Math.round(pur[j].n)+' pairs')+' - '+(pur[j].effect||'NA')});
    }
  }
  return insBox('Selection',narr,chips);
}
function insADNA(){
  var S=(R&&R.samples)||[];
  var anc=[],i;
  for(i=0;i<S.length;i++){ if(S[i].anc) anc.push(S[i]); }
  if(!anc.length){
    return insBox('aDNA','No samples in this run carry ancient-DNA (aDNA) status — nothing to authenticate.');
  }
  var athr=R.anc_thresholds||null;
  var lim=(athr&&athr.damage_min_ct!=null)?athr.damage_min_ct:null;
  var floorTxt=(lim!=null)?(lim*100).toFixed(0)+'%':null;
  var fail=[],nodata=[],pass=[],s,j;
  for(j=0;j<anc.length;j++){
    s=anc[j];
    if(!s.dmg){ nodata.push(s); continue; }
    if(s.f && s.f.indexOf('DAMAGE_LOW')>=0){ fail.push(s); }
    else { pass.push(s); }
  }
  var n=anc.length, plur=(n===1?'':'s'), narr, floorLbl=(floorTxt!=null?'<b>'+floorTxt+'</b> terminal C&gt;T':'terminal C&gt;T authentication');
  if(fail.length){
    narr='<b class="tone-bad">'+fail.length+' of '+n+'</b> aDNA sample'+plur+' claiming antiquity fall below the '+floorLbl+' deamination floor — possible modern DNA mislabelled ancient.';
    if(pass.length) narr+=' <b>'+pass.length+'</b> authenticate.';
    if(nodata.length) narr+=' <b class="tone-warn">'+nodata.length+'</b> lack damage data.';
  } else if(pass.length){
    narr='All <b>'+pass.length+'</b> aDNA sample'+(pass.length===1?'':'s')+' show terminal deamination above the '+(floorTxt!=null?'<b>'+floorTxt+'</b> ':'')+'authentication floor — consistent with post-mortem damage (a screen, not proof of antiquity).';
    if(nodata.length) narr+=' <b class="tone-warn">'+nodata.length+'</b> lack damage data.';
  } else {
    narr='<b class="tone-warn">'+nodata.length+' of '+n+'</b> aDNA sample'+plur+' lack mapDamage profiles — cannot authenticate; run mapDamage2 / pass --mapdamage-dir.';
  }
  var chips=[],a;
  for(a=0;a<fail.length;a++) chips.push({t:fail[a].s,cls:'bad',title:'terminal C>T below the authentication floor'});
  for(a=0;a<nodata.length;a++) chips.push({t:nodata[a].s,cls:'warn',title:'damage data not found'});
  return insBox('aDNA',narr,chips.length?chips:null);
}
function insEpistasis(){
  var E=R.epistasis;
  if(!(E&&E.pairs&&E.pairs.length)) return '';
  var pairs=E.pairs;
  var total=pairs.length;
  var nConc=(E.n_concordant!=null)?E.n_concordant:0;
  var nDisc=(E.n_discordant!=null)?E.n_discordant:0;
  var nStrong=(E.n_strong!=null)?E.n_strong:0;
  // strongest linked pair by |mean r|
  var best=null;
  for(var i=0;i<pairs.length;i++){ var p=pairs[i]; if(best==null||Math.abs(p.r)>Math.abs(best.r)) best=p; }
  function posn(x){ return String(x).split(':').pop(); }
  function lbl(g,pos){ return (g?esc(g):'(intergenic)')+' '+esc(posn(pos)); }
  var arrow=best.direction==='concordant'?'&#8596;':'&#8646;';
  var pairTxt=lbl(best.geneA,best.posA)+' '+arrow+' '+lbl(best.geneB,best.posB);
  var rTxt=(best.r>0?'+':'')+best.r.toFixed(2);
  // headline verdict
  var narr;
  if(nStrong>0){
    narr='<b>'+nStrong+'</b> of '+total+' co-varying variant pair'+(total===1?'':'s')+
      ' clear the FDR bar (<b>strong</b>, q&#8804;0.05) across '+E.n_series+' time series &#8212; '+
      '<b class="tone-warn">'+nConc+' concordant</b> (co-selected / linked) vs <b>'+nDisc+' discordant</b> (competing).';
  } else {
    narr='<b>'+total+'</b> co-varying pair'+(total===1?'':'s')+' surfaced across '+E.n_series+
      ' time series (<b>'+nConc+'</b> concordant, <b>'+nDisc+'</b> discordant), but '+
      '<b class="tone-warn">none reach FDR significance</b> (q&#8804;0.05) &#8212; treat as exploratory.';
  }
  narr+=' Strongest link: <b>'+pairTxt+'</b> at r&#160;=&#160;<b>'+rTxt+'</b>'+
    (best.n>1?(', recurrent in <b>'+best.n+'</b> series'):(' (single series)'))+'.';
  narr+=' Correlation over &#8805;'+E.min_points+' timepoints only flags candidate epistasis &#8212; recurrence across series, not a single trajectory, is what makes it credible.';
  if(E.vars_capped||E.pairs_capped){   // large-cohort scaling caps - say what was covered so nothing is silently dropped
    narr+=' <span class="c">Large cohort: '+
      (E.vars_capped?('scored the top <b>'+E.max_vars+'</b> most-variable variants per series (of up to <b>'+E.n_vars_max+'</b>)'):'')+
      (E.vars_capped&&E.pairs_capped?'; ':'')+
      (E.pairs_capped?('significance-tested the <b>'+E.max_pairs_perm+'</b> strongest of <b>'+E.n_candidates+'</b> candidate pairs'):'')+
      ' &#8212; the full trajectories are in the SNP-dynamics panel / TSV.</span>';
  }
  var chips=[];
  chips.push({t:pairTxt+'  r '+rTxt, cls:best.direction==='concordant'?'warn':'', title:'Strongest co-varying pair &#183; tier '+best.tier+' &#183; perm p '+best.p.toFixed(3)+', FDR q '+best.q.toFixed(3)});
  if(nStrong>0) chips.push({t:nStrong+' strong (q&#8804;0.05)', cls:'warn', title:'Pairs passing Benjamini-Hochberg FDR q &#8804; 0.05'});
  chips.push({t:nConc+' concordant', cls:'', title:'Variants rise and fall together &#8212; candidate linkage / co-selection'});
  chips.push({t:nDisc+' discordant', cls:'', title:'One variant rises as the other falls &#8212; competing lineages / clonal interference'});
  return insBox('Epistasis', narr, chips);
}
function insSnpMatrix(){
  var M=R.snp_matrix;
  if(!(M&&M.rows&&M.rows.length))return '';
  var rows=M.rows, ns=(M.samples&&M.samples.length)?M.samples.length:0;
  if(!ns)return '';
  var nrows=rows.length;
  var total=(M.total_sites!=null)?M.total_sites:nrows;
  var shared=0, singleton=0, maxocc=0, maxrow=null;
  for(var i=0;i<nrows;i++){
    var r=rows[i], occ=(r.n!=null)?r.n:0;   // payload row carries n = occupancy (cells is a sample-index->[af,dp] map, not an array)
    if(occ<=0) continue;
    if(occ===1) singleton++; else shared++;
    if(occ>maxocc){ maxocc=occ; maxrow=r; }
  }
  var classified=shared+singleton;
  if(!classified)return '';
  var sharedPct=Math.round(shared/classified*100);
  var singPct=Math.round(singleton/classified*100);
  var narr='<b>'+fmt(total,'int')+'</b> SNP site'+(total===1?'':'s')+
    ' &#215; <b>'+ns+'</b> sample'+(ns===1?'':'s')+': '+
    '<b>'+fmt(shared,'int')+'</b> ('+sharedPct+'%) shared across &#8805;2 samples, '+
    '<b'+(singPct>=60?' class="tone-warn"':'')+'>'+fmt(singleton,'int')+'</b> ('+singPct+'%) private singletons';
  if(maxrow&&maxocc>=2){
    var gnm=maxrow.gene||maxrow.contig||'';
    narr+=' &#183; largest shared block is <b>'+esc(gnm)+'</b> at '+maxrow.pos+', in <b>'+maxocc+'</b>/'+ns+' samples';
  }
  narr+='.';
  if(singPct>=60) narr+=' Singleton-dominated matrices carry many private calls &#8212; verify low-occupancy sites are real, not per-sample noise.';
  var chips=[];
  chips.push({t:fmt(shared,'int')+' shared', cls:'good', title:'variant sites called in ≥2 samples'});
  chips.push({t:fmt(singleton,'int')+' private', cls:(singPct>=60?'warn':''), title:'variant sites called in exactly 1 sample'});
  if(maxrow&&maxocc>=2){
    chips.push({t:esc(maxrow.gene||maxrow.contig||'?')+' ×'+maxocc, cls:'', title:'most widely shared site: '+esc(maxrow.gene||maxrow.contig||'')+' '+maxrow.pos+' '+esc(maxrow.ref||'')+'→'+esc(maxrow.alt||'')+', in '+maxocc+' of '+ns+' samples'});
  }
  if(M.truncated) chips.push({t:'truncated', cls:'warn', title:'matrix was capped for display ('+fmt(nrows,'int')+' of '+fmt(total,'int')+' sites loaded); shared/private counts are over the loaded sites only — full matrix in the TSV'});
  return insBox('Shared variants', narr, chips);
}
var INS_MAP=[['kraken_body',insKraken],['drug_body',insDrug],['lincomp',insLineages],['plots',insDist],['scatter',insScatter],['corr_body',insCorr],['qcpca_body',insQCspace],['divcomp_body',insRefBias],['genome_body',insGenome],['fn_stacks',insFunction],['gb_body',insGeneBurden],['hot_body',insHotspots],['temporal_body',insTemporal],['pnps_body',insPnps],['adna_body',insADNA],['epi_body',insEpistasis],['snpmx_body',insSnpMatrix]];
