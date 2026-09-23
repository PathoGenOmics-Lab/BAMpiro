// ===================== SUMMARY PAGE: what the run shows, as sentences with the numbers behind them =====================
// Each finding is built from R alone (no DOM), so the same function feeds the summary card, the
// sidebar badge and the panel read-out, and they cannot disagree. A builder returns null when its
// data is absent, and the card simply does not appear.

// The WHO grade 1-2 calls, split into what a whole lineage shares and what only some samples carry.
// A mutation carried by at least DR_WIDE of the samples of a lineage (with at least DR_WIDE_MIN of
// them typed) is that lineage's own: pncA H57D in every M. bovis is why M. bovis resists PZA, and
// counting it as a finding made every M. bovis sample "resistant" and buried the mutations that
// actually arose in the cohort.
var DR_WIDE=0.9, DR_WIDE_MIN=5;
var _drSum=null;
function drSummary(){
  if(_drSum!==null)return _drSum||null;
  var D=R.dr; if(!(D&&D.calls&&D.calls.length)){_drSum=false;return null;}
  var linOf={}; R.samples.forEach(function(s){linOf[s.s]=untyped(s.lineage)?null:linMain(s.lineage);});
  var nLin={}; (D.samples||[]).forEach(function(s){var l=linOf[s];if(l)nLin[l]=(nLin[l]||0)+1;});
  var muts={},order=[];
  D.calls.forEach(function(c){
    if(!(c.gn===1||c.gn===2))return;
    var k=c.drug+'|'+c.gene+'|'+c.mutation, m=muts[k];
    if(!m){m=muts[k]={key:k,drug:c.drug,dr:c.dr,gene:c.gene,mutation:c.mutation,gn:c.gn,grade:c.grade,samples:{}};order.push(k);}
    m.samples[c.s]=(c.af==null?null:c.af);});
  var list=order.map(function(k){var m=muts[k],per={};
    Object.keys(m.samples).forEach(function(s){var l=linOf[s]||'untyped';per[l]=(per[l]||0)+1;});
    m.n=Object.keys(m.samples).length; m.per=per; m.wide=null;
    Object.keys(per).forEach(function(l){if(l!=='untyped'&&(nLin[l]||0)>=DR_WIDE_MIN&&per[l]>=DR_WIDE*nLin[l])m.wide=l;});
    return m;});
  var acquired=list.filter(function(m){return !m.wide;}), wide=list.filter(function(m){return m.wide;});
  var byDrug={},carriers={};
  acquired.forEach(function(m){Object.keys(m.samples).forEach(function(s){carriers[s]=1;(byDrug[m.drug]=byDrug[m.drug]||{})[s]=1;});});
  var drugs=Object.keys(byDrug).map(function(d){return {drug:d,n:Object.keys(byDrug[d]).length,
      muts:acquired.filter(function(m){return m.drug===d;}).sort(function(a,b){return b.n-a.n;})};})
    .sort(function(a,b){return b.n-a.n||(a.drug<b.drug?-1:1);});
  _drSum={all:list,acquired:acquired.sort(function(a,b){return b.n-a.n;}),wide:wide.sort(function(a,b){return b.n-a.n;}),
          carriers:carriers,nCarriers:Object.keys(carriers).length,drugs:drugs,nLin:nLin,nTyped:(D.samples||[]).length};
  return _drSum;
}
function drMutLab(m){return m.gene+' '+m.mutation;}
// "pncA H57D (PZA) in all 82 A4 samples, and in 5 of 93 L7"
function drWideTxt(m,nLin){var others=Object.keys(m.per).filter(function(l){return l!==m.wide&&l!=='untyped';});
  return '<b>'+esc(drMutLab(m))+'</b> ('+esc(m.drug)+') in '+(m.per[m.wide]===nLin[m.wide]?'all ':'')+m.per[m.wide]+(m.per[m.wide]===nLin[m.wide]?'':' of '+nLin[m.wide])+' '+esc(m.wide)+' samples'+
    (others.length?', and in '+others.map(function(l){return m.per[l]+' of '+(nLin[l]||'?')+' '+esc(l);}).join(', '):'');}

// The trajectory screen of the SNP-dynamics panel, over every series at once.
var _dynSum=null;
function dynSummary(){
  if(_dynSum!==null)return _dynSum||null;
  var D=R.dynamics; if(!(D&&D.groups&&D.groups.length)){_dynSum=false;return null;}
  var n=0,sweep=[],rise=0,fall=0,byGene={};
  D.groups.forEach(function(g){g.series.forEach(function(v){n++;var c=dynSelCls(v.traj);
    if(c.dir>0){rise++;var gn=v.gene||'(intergenic)';(byGene[gn]=byGene[gn]||{})[g.group]=1;}else if(c.dir<0)fall++;
    if(c.cls==='sweep')sweep.push({gene:v.gene||'(intergenic)',aa:v.aa,group:g.group,a0:c.a0,aN:c.aN});});});
  var parallel=Object.keys(byGene).map(function(gn){return {gene:gn,n:Object.keys(byGene[gn]).length};})
    .filter(function(x){return x.n>=2;}).sort(function(a,b){return b.n-a.n||(a.gene<b.gene?-1:1);});
  _dynSum={n:n,nSeries:D.groups.length,sweep:sweep,rise:rise,fall:fall,parallel:parallel};
  return _dynSum;
}
// Gene-conversion verdicts, counting one row per event (rep===1) when the cohort step marked them.
// Called events are also split by the QC verdict of the samples carrying them. A tract in a
// sample the QC fails (a mixed culture, a contaminated one) is that sample's problem before it
// is a conversion, and a summary that counts it with the rest overstates what the run found.
// WARN counts with PASS: a sample to review is not a sample to exclude.
function gconvSummary(){
  var G=R.gconv; if(!(G&&G.tracts&&G.tracts.length))return null;
  var rep=G.tracts.some(function(t){return t.rep!=null;}), c={}, samp={};
  var qc={}; R.samples.forEach(function(s){qc[s.s]=s.v;});
  var evPass={}, evFail={}, evBp={};
  G.tracts.forEach(function(t){if(rep&&t.rep!==1)return;c[t.verdict]=(c[t.verdict]||0)+1;
    if(t.verdict==='gene_conversion'){samp[t.s]=1;
      var key=t.event||(t.s+'|'+t.contig+'|'+t.start);
      if(t.bp_reads>0)evBp[key]=1;
      if(qc[t.s]==='FAIL'){(evFail[key]=evFail[key]||{})[t.s]=1;}else evPass[key]=1;}});
  // the failed samples named are those of the events only they carry; a failed sample sharing an
  // event with a kept one is not what "only in samples the QC fails" is about
  var failSamp={};Object.keys(evFail).forEach(function(k){if(!evPass[k])Object.keys(evFail[k]).forEach(function(s){failSamp[s]=1;});});
  var div={}; G.tracts.forEach(function(t){if(t.verdict==='divergent_sample')div[t.s]=1;});
  var nEv=Object.keys(evPass).length+Object.keys(evFail).filter(function(k){return !evPass[k];}).length;
  return {c:c,nCall:c.gene_conversion||0,nSamp:Object.keys(samp).length,bpEvents:Object.keys(evBp).length,nDiv:Object.keys(div).length,rep:rep,
          nEvents:nEv,nPassEvents:Object.keys(evPass).length,failOnly:nEv-Object.keys(evPass).length,failSamples:Object.keys(failSamp).sort()};
}

// ---- the finding cards ----
function _plural(n,one,many){return n===1?one:(many||one+'s');}
// 'Unusually many SNPs' -> 'unusually many SNPs', but an acronym keeps its capitals ('SNPs', 'LoF')
function _lc(t){t=String(t);return (t.length>1&&t.charAt(1)===t.charAt(1).toLowerCase())?t.charAt(0).toLowerCase()+t.slice(1):t;}
function _list(items,max){var a=items.slice(0,max);return a.join(', ')+(items.length>max?' and '+(items.length-max)+' more':'');}
function _tally(list,keep){var t={};list.forEach(function(s){s.f.forEach(function(f){if(keep&&!keep(f))return;t[f]=(t[f]||0)+1;});});
  return Object.keys(t).sort(function(a,b){return t[b]-t[a]||(a<b?-1:1);}).map(function(f){return [f,t[f]];});}
// A FAIL reason as the rule it broke, at the thresholds in force now.
function failPhrase(f){switch(f){
  case 'LOW_DEPTH':return 'depth below '+thr.depth_min+'&#215;';
  case 'LOW_BREADTH':return 'less than '+thr.breadth_min+'% of the genome covered';
  case 'HIGH_MISSING':return 'more than '+thr.missing_max+'% of the consensus missing';
  case 'NO_DATA':return 'no output at all';
  default:return esc(_lc(flagLab(f)));}}
function findQC(){
  var S=R.samples,N=S.length,c=R.counts; if(!N)return null;
  var fail=S.filter(function(s){return s.v==='FAIL';}), warn=S.filter(function(s){return s.v==='WARN';});
  var head='<b>'+c.PASS+'</b> of '+N+' samples are ready to use'+(c.FAIL?'; <b class="tone-bad">'+c.FAIL+'</b> should be excluded':'')+(c.WARN?(c.FAIL?' and ':'; ')+'<b class="tone-warn">'+c.WARN+'</b> need a look':'');
  var body=[];
  if(fail.length)body.push('Excluded for '+_tally(fail,function(f){return FAILF[f];}).map(function(t){return failPhrase(t[0])+' ('+t[1]+')';}).join(', ')+'.');
  if(warn.length)body.push('To review: '+_list(_tally(warn).map(function(t){return esc(_lc(flagLab(t[0])))+' ('+t[1]+')';}),4)+'.');
  if(!fail.length&&!warn.length)body.push('Every sample clears every check at the current thresholds.');
  var w=function(n){return (100*n/N).toFixed(2)+'%';};
  var vis='<div class="f-bar"><span style="width:'+w(c.PASS)+';background:var(--pass)"></span><span style="width:'+w(c.WARN)+';background:var(--warn)"></span><span style="width:'+w(c.FAIL)+';background:var(--fail)"></span></div>';
  return {k:'qc',eyebrow:'Sample QC',tone:c.FAIL?'bad':(c.WARN?'warn':'good'),head:head,body:body.join(' '),vis:vis,
    nums:[{n:c.PASS,l:'ready',t:'good'},{n:c.WARN,l:'to review',t:'warn'},{n:c.FAIL,l:'to exclude',t:'bad'}],
    link:{href:'#flagged',t:'Review the flagged samples'},
    short:c.PASS+' of '+N+' samples are ready'+(c.FAIL?'; '+c.FAIL+' should be excluded':'')+(c.WARN?(c.FAIL?' and ':'; ')+c.WARN+' need a look':'')+'.'};
}
function findIdentity(){
  var K=R.kraken, mism=R.samples.filter(function(s){return s.f.indexOf('LINEAGE_MISMATCH')>=0;});
  var hasK=!!(K&&K.samples&&K.samples.length); if(!hasK&&!mism.length)return null;
  var head='',body=[],nums=[],shortTxt='',tone='good';
  if(hasK){var tgt=krkTarget(),off=[],mix=[];
    K.samples.forEach(function(k){var st0=krkState(k);if(st0==='off')off.push(k);else if(st0==='mixed')mix.push(k);});
    if(off.length){tone='bad';head='<b class="tone-bad">'+off.length+'</b> '+_plural(off.length,'sample is','samples are')+' not <i>'+esc(tgt)+'</i>';
      body.push('Dominated instead by '+_list(krkOffTaxa(off).map(function(t){return t[1]+' &#215; <i>'+esc(t[0])+'</i>';}),4)+'.');
      shortTxt=off.length+' '+_plural(off.length,'sample is','samples are')+' mostly another organism.';}
    else if(mix.length){tone='warn';head='<b class="tone-warn">'+mix.length+'</b> '+_plural(mix.length,'sample carries','samples carry')+' reads from other organisms';
      shortTxt=mix.length+' '+_plural(mix.length,'sample carries','samples carry')+' other organisms.';}
    else{head='Every sample is <i>'+esc(tgt)+'</i>';shortTxt='';}
    if(mix.length)body.push((off.length?mix.length+' more carry':'They carry')+' '+(100-KRK_PURE)+'&#8211;'+(100-KRK_OFF)+'% of other taxa: '+
      _list(mix.map(function(k){return esc(k.s)+(k.secondary?' (<i>'+esc(k.secondary.name)+'</i> '+nv(k.secondary.pct)+'%)':'');}),4)+'.');
    if(!off.length&&!mix.length)body.push('All '+K.samples.length+' samples have at least '+KRK_PURE+'% of their classified reads in it.');
    nums.push({n:K.samples.length-off.length-mix.length,l:'clean',t:'good'},{n:mix.length,l:'mixed',t:'warn'},{n:off.length,l:'other organism',t:'bad'});}
  if(mism.length){tone='bad';
    var t2='<b class="tone-bad">'+mism.length+'</b> '+_plural(mism.length,'sample looks','samples look')+' mapped to the wrong reference';
    var why='They type as a different lineage from the other samples on their reference, so their SNPs measure the distance between the two genomes rather than the samples: '+
      _list(_tally(mism.map(function(s){return {f:[linMain(s.lineage)+' on '+(s.ref||'?')+' ('+linMain(s.ref_lin)+')']};})).map(function(t){return t[1]+' &#215; '+esc(t[0]);}),3)+'.';
    if(head){body.push(t2+'. '+why);}else{head=t2;body.push(why);}
    nums.push({n:mism.length,l:'wrong reference',t:'bad'});
    shortTxt+=(shortTxt?' ':'')+mism.length+' look mapped to the wrong reference.';}
  return {k:'identity',eyebrow:'Identity',tone:tone,head:head,body:body.join(' '),nums:nums,
    link:hasK?{href:'#kraken',t:'See the read composition'}:{href:'#flagged',t:'See the flagged samples'},
    note:hasK?'Kraken2 classifies the reads before mapping; it can only name what its database holds.':'',short:shortTxt};
}
function findLineage(){
  if(!R.lin_present)return null;
  var g={},order=[],N=R.samples.length,notTyped=0;
  R.samples.forEach(function(s){if(untyped(s.lineage)){notTyped++;return;}if(!g[s.lineage]){g[s.lineage]=0;order.push(s.lineage);}g[s.lineage]++;});
  if(!order.length)return null;
  order.sort(function(a,b){return g[b]-g[a];});
  var mixed=R.samples.filter(function(s){return s.f.indexOf('MIXED')>=0;}).length;
  var head=order.length===1?'Every typed sample is <b>'+esc(linLabel(order[0]))+'</b>':
    (order.length===2?'Two lineages: ':order.length+' lineages: ')+_list(order.map(function(l){return '<b>'+esc(linLabel(l))+'</b> ('+g[l]+')';}),3);
  var body=[];
  if(notTyped)body.push(notTyped+' '+_plural(notTyped,'sample','samples')+' could not be typed.');
  body.push(mixed?mixed+' '+_plural(mixed,'sample carries','samples carry')+' a second lineage.':'No sample carries a second lineage.');
  var vis='<div class="f-bar">'+order.map(function(l){return '<span style="width:'+(100*g[l]/N).toFixed(2)+'%;background:'+linColor(l)+'" title="'+esc(linLabel(l))+': '+g[l]+'"></span>';}).join('')+
    (notTyped?'<span style="width:'+(100*notTyped/N).toFixed(2)+'%;background:var(--track)" title="not typed: '+notTyped+'"></span>':'')+'</div>';
  return {k:'lineage',eyebrow:'Lineages',tone:mixed?'warn':'',head:head,body:body.join(' '),vis:vis,link:{href:'#linsum',t:'Lineage summary'}};
}
function findResistance(){
  var d=drSummary(); if(!d)return null;
  var head,body=[],tone='good',shortTxt='';
  if(d.nCarriers){tone='bad';
    head='<b class="tone-bad">'+d.nCarriers+'</b> '+_plural(d.nCarriers,'sample carries','samples carry')+' resistance-associated mutations their lineage does not share';
    body.push(d.drugs.map(function(x){return '<b>'+esc(x.drug)+'</b> in '+x.n+' ('+_list(x.muts.map(function(m){return esc(drMutLab(m))+' '+m.n;}),3)+')';}).join('; ')+'.');
    shortTxt=d.nCarriers+' carry resistance-associated mutations beyond their lineage&#8217;s own ('+d.drugs.map(function(x){return esc(x.drug);}).join(', ')+').';}
  else head='No resistance-associated mutation beyond what each lineage carries';
  if(d.wide.length)body.push('Shared by a whole lineage, so not counted: '+d.wide.slice(0,3).map(function(m){return drWideTxt(m,d.nLin);}).join('; ')+'.');
  return {k:'drug',eyebrow:'Resistance',tone:tone,head:head,body:body.join(' '),
    nums:[{n:d.nCarriers,l:'samples',t:d.nCarriers?'bad':'good'},{n:d.acquired.length,l:'mutations',t:''},{n:d.drugs.length,l:_plural(d.drugs.length,'drug'),t:''}],
    link:{href:'#drug',t:'See every mutation'},note:'WHO catalogue grades 1&#8211;2, called from the reads. A genomic screen, not a drug-susceptibility result.',short:shortTxt};
}
function findDynamics(){
  var d=dynSummary(); if(!d)return null;
  var head='<b>'+d.sweep.length+'</b> '+_plural(d.sweep.length,'allele','alleles')+' swept toward fixation and <b>'+d.rise+'</b> rose, across '+d.nSeries+' series';
  var body=[];
  if(d.parallel.length)body.push('Rising in two or more independent series: '+_list(d.parallel.map(function(x){return '<b>'+esc(x.gene)+'</b> ('+x.n+')';}),5)+'.');
  body.push('Out of '+d.n.toLocaleString('en-US')+' trajectories; '+d.fall+' declined.');
  return {k:'dyn',eyebrow:'Variants over time',tone:'',head:head,body:body.join(' '),
    link:{href:'#dynamics',t:'Open the trajectories'},note:'A screen of allele-frequency trajectories, not a selection test: drift and linked passengers move too.',
    short:d.sweep.length?d.sweep.length+' '+_plural(d.sweep.length,'allele','alleles')+' swept toward fixation.':''};
}
function findGconv(){
  var g=gconvSummary(); if(!g)return null;
  var head='<b>'+g.nPassEvents+'</b> gene-conversion '+_plural(g.nPassEvents,'event')+' called in samples the QC does not fail';
  var body=[];
  if(g.failOnly)body.push(g.failOnly+(g.nPassEvents?' more ':' ')+_plural(g.failOnly,'is','are')+' only in samples the QC fails ('+_list(g.failSamples.map(function(s){return esc(s);}),4)+'), where a mixed or contaminated culture explains donor bases first.');
  // Counted in events like the headline: counted in calls, one event carried by two samples read
  // as more backed events than there were events.
  if(g.bpEvents)body.push((g.bpEvents===g.nEvents?(g.nEvents===1?'That event is':'All '+g.nEvents+' events are'):g.bpEvents+' of the '+g.nEvents+' events '+_plural(g.bpEvents,'is','are'))+' backed by a read crossing a breakpoint.');
  var rest=['ambiguous','coverage_shift','mismapping','reference_artifact'].filter(function(k){return g.c[k];}).map(function(k){return g.c[k]+' '+k.replace(/_/g,' ');});
  if(rest.length)body.push('Also '+rest.join(', ')+'.');
  if(g.nDiv)body.push('<b class="tone-warn">'+g.nDiv+'</b> '+_plural(g.nDiv,'sample was','samples were')+' left out because '+_plural(g.nDiv,'it calls','they call')+' tracts nearly everywhere: '+_plural(g.nDiv,'its genome differs','their genomes differ')+' from the reference as a whole.');
  return {k:'gconv',eyebrow:'Gene conversion',tone:'',head:head,body:body.join(' '),link:{href:'#gconv',t:'Inspect the tracts'},
    note:'Candidates to check in the reads, not confirmed events.'};
}
// Relatedness: the clusters at the pipeline's threshold, and the samples outside their own group.
function relSummary(){var S=R.relatedness;if(!(S&&S.refs))return null;
  var thr=S.threshold,minC=100*(S.min_compared!=null?S.min_compared:0.5),nCl=0,big=0,inCl=0;
  Object.keys(S.refs).forEach(function(k){relClusters(S.refs[k],thr,minC).forEach(function(c){if(c.length>1){nCl++;inCl+=c.length;if(c.length>big)big=c.length;}});});
  return {thr:thr,nCl:nCl,big:big,inCl:inCl,out:R.samples.filter(function(s){return s.grpd;})};}
function findRelatedness(){
  var g=relSummary(); if(!g)return null;
  var head='<b>'+g.nCl+'</b> '+_plural(g.nCl,'cluster')+' of samples within '+g.thr+' SNPs of each other';
  var body=[];
  if(g.nCl)body.push(g.inCl+' samples sit in one, the largest holding '+g.big+'.');
  if(g.out.length)body.push('<b class="tone-warn">'+g.out.length+'</b> '+_plural(g.out.length,'sample sits','samples sit')+' far from the rest of '+_plural(g.out.length,'its group','their groups')+' ('+
    _list(g.out.map(function(s){return esc(s.s);}),4)+'): swapped, mislabelled, contaminated or reinfected.');
  return {k:'rel',eyebrow:'Relatedness',tone:g.out.length?'warn':'',head:head,body:body.join(' '),
    link:{href:'#p-related',t:'See the distances'},
    note:'SNPs between consensus sequences; the threshold is a convention of the organism.',
    short:g.out.length?(g.out.length+' '+_plural(g.out.length,'sample sits','samples sit')+' outside '+_plural(g.out.length,'its group','their groups')+'.'):''};
}
function findCoverage(){
  var S=R.samples; if(!S.length)return null;
  function med(k){return _median(S.map(function(s){return s.m[k];}));}
  var d=med('mean_depth'),b=med('breadth_pct'),c=med('callable_pct'); if(d==null&&b==null)return null;
  var head='Median depth <b>'+(d!=null?fmt(d,'float')+'&#215;':'NA')+'</b>'+(b!=null?', breadth <b>'+b.toFixed(1)+'%</b>':'');
  var rd=_range(S.map(function(s){return s.m.mean_depth;}));
  var body=(rd?'Depth ranges from '+fmt(rd[0],'float')+'&#215; to '+fmt(rd[1],'float')+'&#215;. ':'')+(c!=null?'The median consensus has '+c.toFixed(1)+'% of the genome as confident bases.':'');
  return {k:'cov',eyebrow:'Coverage',tone:'',head:head,body:body,link:{href:'#dist',t:'See the distributions'}};
}
function findings(){return [findQC(),findIdentity(),findRelatedness(),findResistance(),findDynamics(),findLineage(),findGconv(),findCoverage()].filter(function(f){return f;});}
function findingCard(f){
  var nums=(f.nums&&f.nums.length)?'<div class="f-nums">'+f.nums.map(function(x){return '<div class="f-num'+(x.t&&x.n?' '+x.t:'')+'"><b>'+x.n+'</b><span>'+esc(x.l)+'</span></div>';}).join('')+'</div>':'';
  return '<article class="finding'+(f.tone?' f-'+f.tone:'')+'" data-k="'+f.k+'"><div class="f-eyebrow">'+esc(f.eyebrow)+'</div>'+
    '<h3 class="f-head">'+f.head+'</h3>'+(f.vis||'')+(f.body?'<p class="f-body">'+f.body+'</p>':'')+nums+
    '<div class="f-foot">'+(f.note?'<span class="f-note">'+f.note+'</span>':'<span></span>')+
    (f.link?'<a class="f-link" href="'+f.link.href+'">'+esc(f.link.t)+' &#8594;</a>':'')+'</div></article>';
}
function renderExec(){
  var host=el('exec_body'); if(!host)return;
  var F=findings(), S=R.samples;
  var inshort=F.map(function(f){return f.short;}).filter(function(t){return t;}).join(' ');
  var lede=el('sum_lede');
  if(lede)lede.innerHTML=S.length+' samples'+(R.provenance&&R.provenance.reference?' mapped against <b>'+esc(R.provenance.reference)+'</b>':'')+
    ', processed '+esc(String(R.generated||'').split(' ')[0])+'. Every finding links to the page that holds its evidence.';
  host.innerHTML=(inshort?'<div class="inshort"><b>In short.</b> '+inshort+'</div>':'')+
    '<div class="findings">'+F.map(findingCard).join('')+'</div>';
}
// Sidebar badges: the number that says whether a page needs attention.
function renderNavBadges(){
  function set(id,txt,cls,title){var b=el(id);if(!b)return;b.textContent=txt||'';b.className='toc-badge'+(txt?' '+cls:'');b.title=title||'';}
  var c=R.counts;
  set('bdg-qc',c.FAIL?String(c.FAIL):(c.WARN?String(c.WARN):''),c.FAIL?'bad':'warn',c.FAIL?c.FAIL+' samples to exclude':(c.WARN?c.WARN+' samples to review':''));
  var d=drSummary(); set('bdg-drug',d&&d.nCarriers?String(d.nCarriers):'','bad',d&&d.nCarriers?d.nCarriers+' samples with resistance mutations beyond their lineage':'');
  var g=gconvSummary(); set('bdg-gconv',g&&g.nPassEvents?String(g.nPassEvents):'','neu',g?g.nPassEvents+' gene-conversion events in samples the QC does not fail':'');
  var rl=relSummary(); set('bdg-related',rl&&rl.out.length?String(rl.out.length):'','warn',rl&&rl.out.length?rl.out.length+' samples far from the rest of their group':'');
  var y=dynSummary(); set('bdg-variants',y&&y.sweep.length?String(y.sweep.length):'','neu',y?y.sweep.length+' alleles sweeping toward fixation':'');
}
