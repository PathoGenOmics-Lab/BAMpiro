var DYNCOL={fixation:'#2f6fed',emergence:'#1f9d6b',loss:'#e6893a',nonsyn:'#d1495b',high_impact:'#7c3aed'};
var DYNHELP={emergence:'Emergence: the variant is (near-)absent at the first timepoint, then rises above the emergence threshold - a new allele appearing in this series.',fixation:'Fixation: the allele frequency reaches near 1.0 by the last timepoint - the variant has (almost) taken over.',loss:'Loss: the variant is present early then falls back toward 0 - an allele being lost from the series.',nonsyn:'Non-synonymous: the variant changes the protein (missense / stop / frameshift / splice / inframe indel), per snpEff - potentially functional.',high_impact:'High impact: snpEff predicts a HIGH-impact effect (frameshift, stop gained/lost...) - likely to disrupt the gene.'};
var dynState={sel:null,q:'',showAll:false,cap:24,chipsAll:false};
var DYN_CARDS=24, DYN_CHIPS=40;   // cards drawn per step, gene chips listed before "show every gene": 800 cards took 75,000 px
var dynFilter={};
var dynZoom=250;   // trajectory-card width in px (zoom slider); smaller -> more charts per row
var dynShowDP=true;   // draw the per-timepoint read-depth (DP) bars behind each trajectory
function dynHasFlag(v){return v.flags&&v.flags.length;}
function dynColor(flags){ if(!flags)return '#9fb0c3'; if(flags.indexOf('fixation')>=0)return DYNCOL.fixation; if(flags.indexOf('emergence')>=0)return DYNCOL.emergence; if(flags.indexOf('loss')>=0)return DYNCOL.loss; if(flags.indexOf('high_impact')>=0)return DYNCOL.high_impact; return '#5b6b7e'; }
function dynMiniChart(v,th,showDP){
  var n=v.times.length,W=250,H=150,ml=30,mt=10,mb=26;
  var dps=v.dp||[], hasDP=showDP&&dps.some(function(d){return d!=null;});
  var maxDP=1; if(hasDP){ dps.forEach(function(d){ if(d!=null&&d>maxDP)maxDP=d; }); }
  var mr=hasDP?16:10, pw=W-ml-mr, ph=H-mt-mb;
  var bw=Math.min(n<=1?18:(pw/n)*0.5, 16);        // DP bar width (px)
  var px=hasDP?(bw/2+1.5):3, iw=pw-2*px;          // inner x-padding so the first/last point (and its DP bar) clear the value axis
  function X(i){ return ml+px+(n<=1?iw/2:(i/(n-1))*iw); }
  function Y(a){ return mt+(1-a)*ph; }                  // allele frequency (left axis)
  function YD(d){ return mt+ph-(d/maxDP)*ph; }          // read depth (right axis)
  var col=dynColor(v.flags), nonsyn=v.flags&&v.flags.indexOf('nonsyn')>=0;
  var svg='<svg viewBox="0 0 '+W+' '+H+'" width="100%" style="display:block"><title>Allele frequency (0-1, left axis, line) across timepoints'+(hasDP?'; read depth DP as bars with the value on top':'')+'. Hover for exact values.</title>';
  [0,0.5,1].forEach(function(a){ svg+='<line x1="'+ml+'" y1="'+Y(a).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(a).toFixed(1)+'" stroke="'+TH.grid+'" stroke-width="1"/><text x="'+(ml-6)+'" y="'+(Y(a)+3.5).toFixed(1)+'" text-anchor="end" font-size="10" fill="'+TH.mut+'">'+a.toFixed(1)+'</text>'; });
  svg+='<line x1="'+ml+'" y1="'+mt+'" x2="'+ml+'" y2="'+(mt+ph).toFixed(1)+'" stroke="'+TH.axis+'" stroke-width="1"/>';   // left value axis + ticks = a real-figure cue
  [0,0.5,1].forEach(function(a){ svg+='<line x1="'+(ml-3)+'" y1="'+Y(a).toFixed(1)+'" x2="'+ml+'" y2="'+Y(a).toFixed(1)+'" stroke="'+TH.axis+'" stroke-width="1"/>'; });
  if(hasDP){   // depth bars behind the AF line; each bar carries its DP value on top (see the pass after the line)
    dps.forEach(function(d,i){ if(d==null)return; var x=X(i), y=YD(d), h=(mt+ph)-y; svg+='<rect x="'+(x-bw/2).toFixed(1)+'" y="'+y.toFixed(1)+'" width="'+bw.toFixed(1)+'" height="'+Math.max(0,h).toFixed(1)+'" fill="#7ea8d6" opacity="0.45" rx="1.5"><title>t='+esc(v.times[i]==null?i:v.times[i])+'  DP='+d+'</title></rect>'; });
  }
  if(th){ [[th.emerge,DYNCOL.emergence],[th.fix,DYNCOL.fixation]].forEach(function(t){ svg+='<line x1="'+ml+'" y1="'+Y(t[0]).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(t[0]).toFixed(1)+'" stroke="'+t[1]+'" stroke-dasharray="3 3" opacity="0.3"/>'; }); }
  var uid=(window.__qcAF=(window.__qcAF||0)+1), gid='afg'+uid;
  var ptsA=v.traj.map(function(a,i){return X(i).toFixed(1)+','+Y(a).toFixed(1);}), pts=ptsA.join(' '), y0=Y(0).toFixed(1), lastI=v.traj.length-1;
  // vertical gradient wash: event colour strong at the line, dissolving to the baseline (per-chart unique id)
  svg+='<defs><linearGradient id="'+gid+'" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="'+col+'" stop-opacity="0.32"/><stop offset="0.55" stop-color="'+col+'" stop-opacity="0.10"/><stop offset="1" stop-color="'+col+'" stop-opacity="0"/></linearGradient></defs>';
  svg+='<path d="M'+X(0).toFixed(1)+','+y0+' L'+ptsA.join(' L')+' L'+X(lastI).toFixed(1)+','+y0+' Z" fill="url(#'+gid+')"/>';
  svg+='<polyline points="'+pts+'" fill="none" stroke="'+col+'" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>';
  // intermediate points: a panel halo lifts each bead above the line, then the solid dot (nonsyn keeps its ring)
  v.traj.forEach(function(a,i){ if(i===lastI)return; var cx=X(i).toFixed(1), cy=Y(a).toFixed(1); svg+='<circle cx="'+cx+'" cy="'+cy+'" r="3.4" fill="'+TH.panel+'"/><circle cx="'+cx+'" cy="'+cy+'" r="2.4" fill="'+col+'"'+(nonsyn?' stroke="'+DYNCOL.nonsyn+'" stroke-width="1"':'')+'/>'; });
  // endpoint: bold hollow "last value" ring + core + tabular AF label (the finished-figure signature)
  var ex=X(lastI), ey=Y(v.traj[lastI]), ev=v.traj[lastI];
  svg+='<circle cx="'+ex.toFixed(1)+'" cy="'+ey.toFixed(1)+'" r="4.4" fill="'+TH.panel+'" stroke="'+col+'" stroke-width="2.4"/>';
  svg+='<circle cx="'+ex.toFixed(1)+'" cy="'+ey.toFixed(1)+'" r="1.5" fill="'+(nonsyn?DYNCOL.nonsyn:col)+'"/>';
  var labX=ex+7, anchor='start'; if(ex>W-mr-22){ labX=ex-7; anchor='end'; }
  svg+='<text x="'+labX.toFixed(1)+'" y="'+(ey+3.4).toFixed(1)+'" text-anchor="'+anchor+'" font-size="9.5" font-weight="600" fill="'+TH.ink+'" style="font-variant-numeric:tabular-nums">'+ev.toFixed(2)+'</text>';
  // invisible hit targets keep the full per-point AF/DP tooltip on every point despite the thinner beads
  v.traj.forEach(function(a,i){ svg+='<circle cx="'+X(i).toFixed(1)+'" cy="'+Y(a).toFixed(1)+'" r="6" fill="transparent"><title>t='+esc(v.times[i]==null?i:v.times[i])+'  AF='+a.toFixed(3)+(hasDP&&dps[i]!=null?('  DP='+dps[i]):'')+'</title></circle>'; });
  if(hasDP){ dps.forEach(function(d,i){ if(d==null)return; svg+='<text x="'+X(i).toFixed(1)+'" y="'+(YD(d)-3).toFixed(1)+'" text-anchor="middle" font-size="8.5" font-weight="600" fill="'+TH.ink+'" stroke="'+TH.panel+'" stroke-width="2.6" paint-order="stroke" style="paint-order:stroke">'+d+'</text>'; }); }   // DP value on top of each bar
  v.times.forEach(function(t,i){ svg+='<text x="'+X(i).toFixed(1)+'" y="'+(H-8)+'" text-anchor="middle" font-size="10.5" fill="'+TH.mut+'">'+esc(t==null?i:t)+'</text>'; });
  svg+='</svg>';
  return svg;
}
// ---- Selection screen: turn the AF trajectories into an analytical read-out ----
// For each trajectory we fit logit(AF) vs time (OLS) to get an apparent selection coefficient s,
// classify the move (sweep / emerging / rising / declining / lost / stable), cross-reference the
// drug-resistance catalogue by gene+mutation, and flag genes rising in >=2 independent series
// (candidate convergent/parallel adaptation). It is a heuristic screen, NOT a formal selection test.
function dynLogit(p){ p=Math.max(0.02,Math.min(0.98,p)); return Math.log(p/(1-p)); }
function dynSlope(traj,times){ var xs=[],ys=[],i; for(i=0;i<traj.length;i++){ if(traj[i]!=null){ var t=(times&&times[i]!=null&&!isNaN(+times[i]))?+times[i]:i; xs.push(t); ys.push(dynLogit(traj[i])); } }
  var n=xs.length; if(n<2) return null; var mx=0,my=0; for(i=0;i<n;i++){mx+=xs[i];my+=ys[i];} mx/=n;my/=n;
  var sxy=0,sxx=0,syy=0; for(i=0;i<n;i++){ var dx=xs[i]-mx,dy=ys[i]-my; sxy+=dx*dy; sxx+=dx*dx; syy+=dy*dy; }
  if(sxx<=0) return {s:0,r2:0,n:n}; return {s:sxy/sxx, r2:(syy>0?(sxy*sxy)/(sxx*syy):1), n:n}; }
function dynSelCls(traj){ var vals=[],i; for(i=0;i<traj.length;i++){ if(traj[i]!=null) vals.push(traj[i]); }
  if(vals.length<2) return {cls:'single',dir:0,delta:0,a0:vals[0]||0,aN:vals[0]||0};
  var a0=vals[0],aN=vals[vals.length-1],delta=aN-a0;
  var cls='stable',dir=0;
  if(delta>=0.15&&a0<=0.25&&aN>=0.75){cls='sweep';dir=1;}
  else if(delta>=0.15&&a0<=0.1){cls='emerge';dir=1;}
  else if(delta>=0.15){cls='rising';dir=1;}
  else if(delta<=-0.15&&aN<=0.15){cls='lost';dir=-1;}
  else if(delta<=-0.15){cls='declining';dir=-1;}
  return {cls:cls,dir:dir,delta:delta,a0:a0,aN:aN}; }
var DYNCLS={sweep:{lab:'sweep → fixation',col:'#2f6fed'},emerge:{lab:'emerging',col:'#1f9d6b'},rising:{lab:'rising',col:'#2ea36b'},declining:{lab:'declining',col:'#e6893a'},lost:{lab:'lost',col:'#e0544f'},stable:{lab:'stable',col:'#8895a6'},single:{lab:'single point',col:'#8895a6'}};
function dynDRindex(){ var idx={}; if(R.dr&&R.dr.calls){ R.dr.calls.forEach(function(c){ if(!c.gene)return; var k=((c.gene||'')+'|'+(c.mutation||'')).toLowerCase().replace(/\s+/g,''); if(!idx[k]||((c.gn===1||c.gn===2)&&!(idx[k].gn===1||idx[k].gn===2))) idx[k]=c; }); } return idx; }
function dynDRmatch(idx,v){ if(!v.aa) return null; var k=((v.gene||'')+'|'+v.aa).toLowerCase().replace(/\s+/g,''); return idx[k]||null; }
function dynSpark(traj,col){ var n=traj.length; if(n<2) return ''; var W=110,H=28,pad=3;
  function X(i){return pad+(i/(n-1))*(W-2*pad);} function Y(a){return H-pad-a*(H-2*pad);}
  var pts=[],i; for(i=0;i<n;i++){ if(traj[i]!=null) pts.push(X(i).toFixed(1)+','+Y(traj[i]).toFixed(1)); }
  if(pts.length<2) return ''; var y0=Y(0).toFixed(1), li=n-1;
  return '<svg width="'+W+'" height="'+H+'" style="display:block"><path d="M'+pts.join(' L')+' L'+X(li).toFixed(1)+','+y0+' L'+X(0).toFixed(1)+','+y0+' Z" fill="'+col+'" opacity="0.15"/>'+
    '<polyline points="'+pts.join(' ')+'" fill="none" stroke="'+col+'" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"/>'+
    '<circle cx="'+X(li).toFixed(1)+'" cy="'+Y(traj[li]).toFixed(1)+'" r="2.5" fill="'+col+'"/></svg>'; }
function renderDynamics(){
  var host=el('dyn_body'), sec=el('dynamics'); if(!host)return;
  var D=R.dynamics;
  if(!(D&&D.groups&&D.groups.length)){ if(sec)sec.style.display='none'; return; }
  if(sec)sec.style.display='';
  var th=D.thresholds||{emerge:0.25,fix:0.9,loss:0.1};
  var vars=[];
  D.groups.forEach(function(g){ g.series.forEach(function(s){ vars.push({gene:s.gene||'(intergenic)',pos:s.pos,group:g.group,times:g.times,traj:s.traj,dp:s.dp,flags:s.flags,eff:s.eff,imp:s.imp,alt:s.alt,aa:s.aa,aa_h37rv:s.aa_h37rv,pos_h37rv:s.pos_h37rv}); }); });
  // per-series (group) metadata + the fields usable as a series filter: group-invariant (one value per
  // series, so timepoint/date -> the trajectory axis -> excluded) and with >1 value across series.
  var groupMeta={}; D.groups.forEach(function(g){ groupMeta[g.group]=g.meta||{}; });
  var mfields=(R.sample_meta&&R.sample_meta.fields)||[];
  var dynFields=[], dynFieldVals={};
  mfields.forEach(function(f){
    var invariant=true, all={};
    D.groups.forEach(function(g){ var vals=(g.meta&&g.meta[f])||[]; if(vals.length>1) invariant=false; vals.forEach(function(v){all[v]=1;}); });
    var distinct=Object.keys(all);
    if(invariant && distinct.length>1){ dynFields.push(f); dynFieldVals[f]=distinct.sort(); }
  });
  function passFilter(grp){ var m=groupMeta[grp]||{}; for(var f in dynFilter){ if(dynFilter[f]){ var vals=m[f]||[]; if(vals.indexOf(dynFilter[f])<0) return false; } } return true; }
  function filterActive(){ for(var f in dynFilter){ if(dynFilter[f]) return true; } return false; }
  var genes, geneList, singleGroup=false, visGroupName='';
  function recompute(){
    genes={}; var gset={};
    vars.forEach(function(v){ if(passFilter(v.group)){ (genes[v.gene]=genes[v.gene]||[]).push(v); gset[v.group]=1; } });
    var gkeys=Object.keys(gset);
    singleGroup=gkeys.length<=1;   // one visible series -> the per-card group tag is redundant
    visGroupName=gkeys.length===1?gkeys[0]:'';
    geneList=Object.keys(genes).sort(function(a,b){
      var fa=genes[a].filter(dynHasFlag).length, fb=genes[b].filter(dynHasFlag).length;
      return fb-fa || genes[b].length-genes[a].length || a.localeCompare(b);
    });
  }
  recompute();
  if(dynState.sel===null){
    dynState.sel={};
    var fg=geneList.filter(function(g){return genes[g].some(dynHasFlag);});
    (fg.length?fg:geneList).slice(0,8).forEach(function(g){dynState.sel[g]=1;});   // collapsed default: the top few genes; 'show all' reveals the rest
  }
  var filterUI=dynFields.length?('<div class="dyn-filters"><span class="snpmx-flabel" title="Show only the connected series (patient / passage line...) matching these metadata values. The time axis of each trajectory is unchanged.">filter series:</span>'+
    dynFields.map(function(f){ return '<label class="snpmx-fsel">'+esc(f)+' <select data-df="'+esc(f)+'"><option value="">all</option>'+
      dynFieldVals[f].map(function(v){return '<option value="'+esc(v)+'"'+(dynFilter[f]===v?' selected':'')+'>'+esc(v)+'</option>';}).join('')+'</select></label>'; }).join('')+
    '<button class="dyn-btn" id="dynfclear">clear</button></div>'):'';
  host.innerHTML=
    '<div class="insight dyn-insight" id="dyn_insight"></div>'+
    '<div class="dyn-controls">'+
      '<input id="dynsearch" class="dyn-search" type="search" title="Type a gene name to filter the gene chips and the grid below" placeholder="search gene..." value="'+esc(dynState.q)+'">'+
      '<button class="dyn-btn" id="dynFlag" title="Show only genes that have at least one flagged variant">flagged genes</button>'+
      '<button class="dyn-btn showall-btn" id="dynAll" title="Show every moving SNP&#39;s trajectory (select all genes)">show all</button>'+
      '<button class="dyn-btn" id="dynNone" title="Deselect all genes">clear</button>'+
      '<label class="dyn-toggle" title="Show a per-timepoint read-depth (DP) bar behind each trajectory"><input type="checkbox" id="dynDP"'+(dynShowDP?' checked':'')+'> depth bars</label>'+
      '<label class="dyn-zoom" title="Resize the trajectory cards - drag left to fit more charts per row"><span>'+icon('search','sort')+'&#8211;/+</span><input type="range" id="dynzoom" min="165" max="360" step="5" value="'+dynZoom+'"></label>'+
      '<span class="dyn-count" id="dynCount"></span></div>'+
    filterUI+
    '<div class="dyn-legend"><span title="'+DYNHELP.emergence+'"><i style="background:'+DYNCOL.emergence+'"></i>emergence <span class="infoi">i</span></span><span title="'+DYNHELP.fixation+'"><i style="background:'+DYNCOL.fixation+'"></i>fixation <span class="infoi">i</span></span><span title="'+DYNHELP.loss+'"><i style="background:'+DYNCOL.loss+'"></i>loss <span class="infoi">i</span></span><span title="'+DYNHELP.nonsyn+'"><i style="border:2px solid '+DYNCOL.nonsyn+';background:var(--panel)"></i>non-synonymous <span class="infoi">i</span></span></div>'+
    '<div class="dyn-genechips" id="dynchips"></div>'+
    '<div class="dyn-grid" id="dyngrid"></div>';
  function paintChips(){
    var q=dynState.q.toLowerCase();
    var all=geneList.filter(function(g){return !q||g.toLowerCase().indexOf(q)>=0;});
    var capped=!dynState.chipsAll&&all.length>DYN_CHIPS, list=capped?all.filter(function(g,i){return i<DYN_CHIPS||dynState.sel[g];}):all;
    el('dynchips').innerHTML=list.length?list.map(function(g){
      var nf=genes[g].filter(dynHasFlag).length;
      return '<button class="dyn-chip'+(dynState.sel[g]?' sel':'')+'" data-g="'+esc(g)+'" title="Click to toggle. '+genes[g].length+' variant trajectory(ies)'+(nf?(', '+nf+' flagged'):'')+'">'+esc(g)+' <b>'+genes[g].length+'</b>'+(nf?'<i class="dyn-dot" title="'+nf+' flagged variant(s) in this gene"></i>':'')+'</button>';
    }).join('')+(capped?'<button class="dyn-btn" id="dynChipsAll" title="List every gene with a trajectory">+'+(all.length-list.length)+' more genes</button>':''):'<span class="c">'+(filterActive()?'no gene matches the current series filter':'no gene matches "'+esc(dynState.q)+'"')+'</span>';
    Array.prototype.forEach.call(el('dynchips').querySelectorAll('.dyn-chip'),function(b){ b.onclick=function(){ var g=b.getAttribute('data-g'); if(dynState.sel[g])delete dynState.sel[g]; else dynState.sel[g]=1; dynState.cap=DYN_CARDS; paintChips(); paintGrid(); }; });
    var ca=el('dynChipsAll'); if(ca)ca.onclick=function(){ dynState.chipsAll=true; paintChips(); };
  }
  function paintGrid(){
    var q=dynState.q.toLowerCase();
    var sel=geneList.filter(function(g){return dynState.sel[g] && (!q||g.toLowerCase().indexOf(q)>=0);});
    var nvar=0; sel.forEach(function(g){nvar+=genes[g].length;});
    el('dynCount').innerHTML=sel.length+' of '+geneList.length+' genes'+(sel.length?(' &#183; '+nvar+' trajectories'):'')+(singleGroup&&visGroupName?(' &#183; series <b>'+esc(visGroupName)+'</b>'):'')+(filterActive()?' (filtered)':'');
    var grid=el('dyngrid');
    grid.style.setProperty('--dyncw', dynZoom+'px');
    if(!sel.length){ grid.innerHTML='<div class="dyn-empty" style="grid-column:1/-1">&#128204; '+(geneList.length?(dynState.q?('no selected gene matches &quot;'+esc(dynState.q)+'&quot;'):'Search and select one or more genes above to see the allele-frequency trajectories of their variants.'):'no variant trajectory matches the current series filter.')+'</div>'; return; }
    var cards=[], total=nvar;
    sel.forEach(function(g){
      genes[g].slice().sort(function(a,b){ return ((dynHasFlag(b)?1:0)-(dynHasFlag(a)?1:0)) || (String(a.pos)>String(b.pos)?1:-1); }).forEach(function(v){
        var flagged=dynHasFlag(v);
        var chips=(v.flags||[]).map(function(f){return '<span class="dyn-fchip" style="background:'+(DYNCOL[f]||'#8895a6')+'" title="'+(DYNHELP[f]||f)+'">'+f+'</span>';}).join('');
        var posNum=String(v.pos).split(':').pop();
        cards.push('<div class="dyn-card'+(flagged?' flagged':'')+'">'+
          '<div class="dyn-card-h"><span class="dyn-cardgene" title="Gene (click its chip above to toggle)">'+esc(v.gene||'(intergenic)')+'</span>'+geneRvTag(v.gene)+(singleGroup?'':'<span class="dyn-grp" title="Connected series this variant belongs to (the metadata group column, e.g. patient / passage line)">'+esc(v.group)+'</span>')+'<span class="dyn-pos" title="Genomic position (contig:position) of this SNP: '+esc(v.pos)+(v.pos_h37rv?('  &#183; '+esc(AA2LBL)+': '+esc(v.pos_h37rv)):'')+'">'+refPos(posNum,v.pos_h37rv)+'</span></div>'+
          '<div class="dyn-eff" title="Predicted effect (snpEff) and protein change HGVS.p: ref amino acid, codon position, alt amino acid">'+esc(v.eff||'variant')+(v.aa?(' &#183; <b class="dyn-aa">'+aaDual(v.aa,v.aa_h37rv)+'</b>'):(v.alt?(' &#183; &#8594;'+esc(v.alt)):''))+'</div>'+
          dynMiniChart(v,th,dynShowDP)+
          '<div class="dyn-card-f">'+(chips||'<span class="c" title="no emergence / fixation / loss / non-synonymous event for this variant">no event</span>')+'</div>'+
        '</div>');
      });
    });
    var shown=Math.min(cards.length,dynState.cap);
    grid.innerHTML=cards.slice(0,shown).join('')+(cards.length>shown?'<div class="dyn-more" style="grid-column:1/-1"><button class="dyn-btn" id="dynMore">show '+Math.min(DYN_CARDS,cards.length-shown)+' more</button><span class="c">'+shown+' of '+total+' trajectories drawn; hover a point for its allele frequency and depth</span></div>':'');
    var mb=el('dynMore'); if(mb)mb.onclick=function(){ dynState.cap+=DYN_CARDS; paintGrid(); };
  }
  function paintInsight(){
    var box=el('dyn_insight'); if(!box) return;
    var drIdx=dynDRindex();
    var scored=vars.filter(function(v){ return passFilter(v.group); }).map(function(v){ var sl=dynSlope(v.traj,v.times), c=dynSelCls(v.traj); return {v:v,s:sl?sl.s:0,r2:sl?sl.r2:0,cls:c.cls,dir:c.dir,delta:c.delta,a0:c.a0,aN:c.aN,dr:dynDRmatch(drIdx,v)}; });
    // what rises comes first (sweeps, then the largest gains), then what falls: under a drug the gains are the story
    var CLSRANK={sweep:0,emerge:1,rising:2,declining:3,lost:4};
    var movers=scored.filter(function(x){ return x.dir!==0; }).sort(function(a,b){ return ((CLSRANK[a.cls]<3?0:1)-(CLSRANK[b.cls]<3?0:1))||(Math.abs(b.delta)-Math.abs(a.delta))||(Math.abs(b.s)-Math.abs(a.s)); });
    if(!movers.length){ box.innerHTML='<div class="dyn-ins-h"><span class="dyn-ins-title">Selection screen</span></div><div class="dyn-ins-narr">No trajectory shows a directional allele-frequency change beyond noise'+(filterActive()?' in the current series filter':'')+' &#8212; the alleles present look static across the sampled timepoints.</div>'; return; }
    var sweeps=movers.filter(function(x){return x.cls==='sweep';});
    var drUp=movers.filter(function(x){return x.dr&&x.dir>0;});
    var drugsUp={}; drUp.forEach(function(x){ if(x.dr.drug) drugsUp[x.dr.drug]=1; });
    var byGene={}; movers.filter(function(x){return x.dir>0;}).forEach(function(x){ var g=x.v.gene||'(intergenic)'; (byGene[g]=byGene[g]||{})[x.v.group]=1; });
    var conv=[]; for(var g in byGene){ var ser=Object.keys(byGene[g]); if(ser.length>=2) conv.push({gene:g,series:ser.sort()}); }
    conv.sort(function(a,b){return b.series.length-a.series.length||(a.gene<b.gene?-1:1);});
    var risers=movers.filter(function(x){return x.dir>0;}).length, fallers=movers.length-risers;
    var narr='<b>'+movers.length+'</b> of <b>'+scored.length+'</b> trajectories are moving directionally &#8212; <b>'+risers+'</b> rising'+(fallers?', <b>'+fallers+'</b> declining':'')+
      (sweeps.length?', <b>'+sweeps.length+'</b> sweeping toward fixation':'')+
      (drUp.length?'. <b class="tone-bad">'+drUp.length+'</b> rising allele'+(drUp.length>1?'s are':' is a')+' known resistance mutation'+(drUp.length>1?'s':'')+' ('+esc(Object.keys(drugsUp).join(', '))+')':'')+
      (conv.length?'. <b class="tone-warn">'+conv.length+'</b> gene'+(conv.length>1?'s':'')+' rising in parallel across independent series &#8212; candidate convergent adaptation':'')+'.';
    var CONVMAX=12, convHTML=conv.length?'<div class="dyn-conv">'+conv.slice(0,CONVMAX).map(function(c){ return '<span class="dyn-conv-chip" title="'+esc(c.gene)+' has a rising variant in '+c.series.length+' independent series ('+esc(c.series.join(', '))+'). The same gene under selection in parallel is a strong signal of real adaptation (e.g. drug pressure), not noise.">'+esc(c.gene)+' &#8593; '+c.series.length+' series</span>'; }).join('')+
      (conv.length>CONVMAX?'<span class="dyn-conv-more">+'+(conv.length-CONVMAX)+' more genes rising in 2 or more series</span>':'')+'</div>':'';
    var rows=movers.slice(0,8).map(function(x){ var v=x.v, cc=DYNCLS[x.cls]||DYNCLS.stable, spCol=(x.dir>0?cc.col:'#e6893a');
      var drTag=x.dr?('<span class="dyn-drtag'+((x.dr.gn===1||x.dr.gn===2)?' r':'')+'" title="Drug-resistance catalogue match: '+esc(x.dr.drug||'')+' &#183; grade '+esc(x.dr.grade||'')+'">'+esc(x.dr.drug||'DR')+'</span>'):'';
      return '<tr class="dyn-ins-row" data-g="'+esc(v.gene||'')+'" title="Click to show this gene in the charts below &#183; logit-slope s='+x.s.toFixed(3)+', fit R&#178;='+x.r2.toFixed(2)+'">'+
        '<td class="dyn-ins-v"><b>'+esc(v.gene||'(intergenic)')+'</b>'+(v.aa?' <span class="dyn-ins-mut">'+esc(v.aa)+'</span>':'')+'<span class="dyn-ins-grp">'+esc(v.group)+'</span></td>'+
        '<td class="dyn-ins-sp">'+dynSpark(v.traj,spCol)+'</td>'+
        '<td class="dyn-ins-tr"><span class="dyn-ins-arrow" style="color:'+spCol+'">'+(x.dir>0?'&#8593;':'&#8595;')+'</span> '+(x.delta>0?'+':'')+x.delta.toFixed(2)+'<span class="dyn-ins-ep">'+x.a0.toFixed(2)+'&#8594;'+x.aN.toFixed(2)+'</span></td>'+
        '<td><span class="dyn-ins-chip" style="background:'+cc.col+'">'+cc.lab+'</span>'+drTag+'</td>'+
      '</tr>';
    }).join('');
    box.innerHTML='<div class="dyn-ins-h"><span class="dyn-ins-title">Selection screen</span><span class="dyn-ins-sub">which alleles are changing, how fast, and whether it looks like selection</span></div>'+
      '<div class="dyn-ins-narr">'+narr+'</div>'+convHTML+
      '<div class="dyn-ins-scroll"><table class="dyn-ins-tbl"><thead><tr><th>variant &#183; series</th><th>trajectory</th><th>change</th><th>call</th></tr></thead><tbody>'+rows+'</tbody></table></div>'+
      '<div class="dyn-ins-caveat">Heuristic screen from a logit-AF slope, not a formal selection test: few timepoints, allele frequencies carry depth noise, and drift or hitchhiking (linkage) can mimic selection. Convergence across independent series is the most robust signal.</div>';
    Array.prototype.forEach.call(box.querySelectorAll('.dyn-ins-row'),function(r){ r.onclick=function(){ var g=r.getAttribute('data-g'); if(!g)return; dynState.sel={}; dynState.sel[g]=1; dynState.q=''; dynState.cap=DYN_CARDS; el('dynsearch').value=''; paintChips(); paintGrid(); el('dyngrid').scrollIntoView({behavior:'smooth',block:'nearest'}); }; });
  }
  el('dynsearch').oninput=function(){ dynState.q=this.value; dynState.cap=DYN_CARDS; paintChips(); paintGrid(); };
  el('dynFlag').onclick=function(){ dynState.sel={}; dynState.cap=DYN_CARDS; geneList.filter(function(g){return genes[g].some(dynHasFlag);}).forEach(function(g){dynState.sel[g]=1;}); paintChips(); paintGrid(); };
  function dynSetAll(on){ dynState.showAll=on; dynState.sel={}; dynState.cap=DYN_CARDS;
    if(on){ geneList.forEach(function(g){dynState.sel[g]=1;}); }
    else { var fg=geneList.filter(function(g){return genes[g].some(dynHasFlag);}); (fg.length?fg:geneList).slice(0,8).forEach(function(g){dynState.sel[g]=1;}); }
    var bb=el('dynAll'); if(bb){ bb.textContent=on?'show less':'show all'; bb.classList.toggle('on',on); }
    paintChips(); paintGrid(); if(window.__syncSitesBtn)window.__syncSitesBtn(); }
  el('dynAll').onclick=function(){ dynSetAll(!dynState.showAll); };   // toggle: every trajectory <-> flagged genes only
  window.__dynSetAll=dynSetAll;
  if(dynState.showAll){ el('dynAll').textContent='show less'; el('dynAll').classList.add('on'); }
  el('dynNone').onclick=function(){ dynState.sel={}; paintChips(); paintGrid(); };
  el('dynzoom').oninput=function(){ dynZoom=+this.value; el('dyngrid').style.setProperty('--dyncw', dynZoom+'px'); };
  el('dynDP').onchange=function(){ dynShowDP=this.checked; paintGrid(); };
  if(dynFields.length){
    Array.prototype.forEach.call(host.querySelectorAll('.dyn-filters select'),function(sel){ sel.onchange=function(){ var f=sel.getAttribute('data-df'); if(sel.value)dynFilter[f]=sel.value; else delete dynFilter[f]; recompute(); paintInsight(); paintChips(); paintGrid(); }; });
    el('dynfclear').onclick=function(){ dynFilter={}; Array.prototype.forEach.call(host.querySelectorAll('.dyn-filters select'),function(s){s.value='';}); recompute(); paintInsight(); paintChips(); paintGrid(); };
  }
  paintInsight(); paintChips(); paintGrid();
}


