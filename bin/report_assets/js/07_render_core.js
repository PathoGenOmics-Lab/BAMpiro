function renderInsights(){
  for(var i=0;i<INS_MAP.length;i++){
    var bodyId=INS_MAP[i][0], fn=INS_MAP[i][1], body=el(bodyId); if(!body||!fn)continue;
    var sec=body.closest?body.closest('section'):null; if(!sec)continue;
    var slotId='ins_'+bodyId, slot=el(slotId);
    if(!slot){ slot=document.createElement('div'); slot.id=slotId; var h2=sec.querySelector('h2'); if(h2&&h2.parentNode){ h2.parentNode.insertBefore(slot,h2.nextSibling); } else { sec.insertBefore(slot,sec.firstChild); } }
    var html=''; try{ html=fn(); }catch(e){ html=''; }
    slot.innerHTML=html||'';
  }
}
function renderExec(){
  var host=el('exec_body'); if(!host)return;
  var S=R.samples, N=S.length, c=R.counts, toEx=c.FAIL;
  function med(k){return _median(S.map(function(s){return s.m[k];}));}
  function rg(k){return _range(S.map(function(s){return s.m[k];}));}
  var mDepth=med('mean_depth'), mBreadth=med('breadth_pct'), rDepth=rg('mean_depth'), rBreadth=rg('breadth_pct');
  var contam=(R.kraken&&R.kraken.samples)?R.kraken.samples.filter(function(k){return k.primary&&k.primary.pct<90;}).length:null;
  var resSamp=null, resDrugs=[];
  if(R.dr&&R.dr.calls){var rs={},dd={}; R.dr.calls.forEach(function(cl){if(cl.gn===1||cl.gn===2){rs[cl.s]=1; if(cl.drug)dd[cl.drug]=(dd[cl.drug]||0)+1;}});
    resSamp=Object.keys(rs).length; resDrugs=Object.keys(dd).sort(function(a,b){return dd[b]-dd[a];}).slice(0,4);}
  function card(l,n,s,tone){return '<div class="kpi'+(tone?' '+tone:'')+'"><div class="kpi-l">'+l+'</div><div class="kpi-n">'+n+'</div><div class="kpi-s">'+(s||'')+'</div></div>';}
  var passRate=N?Math.round(c.PASS/N*100):0;
  var cards=[
    card('Samples', N, '<span class="v PASS xs">'+c.PASS+' pass</span> <span class="v WARN xs">'+c.WARN+' warn</span> <span class="v FAIL xs">'+c.FAIL+' fail</span>'),
    card('Pass rate', passRate+'%', toEx?('<b>'+toEx+'</b> to exclude'):'all usable', passRate>=80?'good':(passRate>=50?'warn':'bad')),
    card('Median depth', (mDepth!=null?fmt(mDepth,'float')+'&#215;':'NA'), rDepth?(fmt(rDepth[0],'float')+'-'+fmt(rDepth[1],'float')+'&#215; range'):''),
    card('Median breadth', (mBreadth!=null?mBreadth.toFixed(1)+'%':'NA'), rBreadth?(rBreadth[0].toFixed(0)+'-'+rBreadth[1].toFixed(0)+'% range'):'')];
  if(contam!=null) cards.push(card('Contamination', contam, contam?'sample(s) &lt; 90% primary':'none flagged', contam?'warn':'good'));
  if(resSamp!=null) cards.push(card('Resistance', resSamp, resDrugs.length?esc(resDrugs.join(' · ')):'no R calls', resSamp?'warn':''));
  var flagc={}; S.forEach(function(s){(s.f||[]).forEach(function(f){flagc[f]=(flagc[f]||0)+1;});});
  var flags=Object.keys(flagc).sort(function(a,b){return flagc[b]-flagc[a];});
  var prof=flags.length?flags.map(function(f){var n=flagc[f],w=Math.round(n/N*100),fatal=FAILF[f];
    return '<div class="qcp-row"><span class="qcp-lab">'+f+'</span><span class="qcp-bar"><span style="width:'+Math.max(4,w)+'%;background:'+(fatal?'var(--fail)':'var(--warn)')+'"></span></span><span class="qcp-n">'+n+'</span></div>';}).join(''):'<div class="krk-mut">No sample trips any check at the current thresholds.</div>';
  var linc={}; S.forEach(function(s){if(s.lineage)linc[s.lineage]=(linc[s.lineage]||0)+1;});
  var lins=Object.keys(linc).sort();
  var strip=lins.length?'<div class="lincomp-bar" style="margin:0 0 8px">'+lins.map(function(l){return '<div class="lseg" style="width:'+(linc[l]/N*100).toFixed(2)+'%;background:'+linColor(l)+'" title="'+esc(l)+' n='+linc[l]+'"></div>';}).join('')+'</div><div class="lincomp-lab" style="padding:0">'+lins.map(function(l){return '<span class="lchip"><i style="background:'+linColor(l)+'"></i>'+esc(l)+' <b>'+linc[l]+'</b></span>';}).join('')+'</div>':'<div class="krk-mut">no lineage calls</div>';
  var narr='<b>'+N+'</b> samples analysed against '+(R.provenance&&R.provenance.reference?esc(R.provenance.reference):'the reference')+' &#183; <b>'+c.PASS+'</b> pass, '+(toEx?'<b class="tone-bad">'+toEx+'</b> recommended for exclusion':'0 to exclude')+' &#183; median depth <b>'+(mDepth!=null?fmt(mDepth,'float')+'&#215;':'NA')+'</b>, breadth <b>'+(mBreadth!=null?mBreadth.toFixed(1)+'%':'NA')+'</b>'+(contam?' &#183; <b class="tone-warn">'+contam+'</b> possibly contaminated':'')+(resSamp?' &#183; drug resistance in <b class="tone-warn">'+resSamp+'</b> sample(s)':'')+'.';
  host.innerHTML='<div class="exec-narr">'+narr+'</div><div class="exec-grid">'+cards.join('')+'</div>'+
    '<div class="exec-cols"><div class="exec-block"><div class="exec-h">QC quality profile <span class="krk-mut">- samples tripping each check (red = gate-failing)</span></div>'+prof+'</div>'+
    '<div class="exec-block"><div class="exec-h">Cohort lineages</div>'+strip+'</div></div>';
}
function renderOverview(){
  var c=R.counts;
  el('summary').innerHTML=donut(c)+'<div class="counts">'+
    '<div class="c all"><div class="n">'+R.samples.length+'</div><div class="l">samples</div></div>'+
    ['PASS','WARN','FAIL'].map(function(v){return '<div class="c '+v.toLowerCase()+'"><div class="n">'+c[v]+'</div><div class="l">'+v.toLowerCase()+'</div></div>';}).join('')+'</div>';
  var freq={}; R.samples.forEach(function(s){s.f.forEach(function(f){freq[f]=(freq[f]||0)+1;});});
  var keys=Object.keys(freq).sort(function(a,b){return freq[b]-freq[a];});
  el('chips').innerHTML='<span class="t">flags</span>'+(keys.length?keys.map(function(f){return '<span class="chip'+(st.flagFilter==f?' on':'')+'" data-f="'+f+'" role="button" tabindex="0" aria-pressed="'+(st.flagFilter==f?'true':'false')+'" aria-label="filter by '+f+'">'+f+'<span class="k">'+freq[f]+'</span></span>';}).join('')+'<span class="chip-hint">click to filter</span>':'<span style="color:#94a3b8;font-size:12px">none - every sample clear ✓</span>');
  Array.prototype.forEach.call(document.querySelectorAll('#chips .chip'),function(ch){function tog(){var f=ch.getAttribute('data-f');st.flagFilter=(st.flagFilter==f?null:f);st.onlyFlagged=false;renderAll();}
    ch.onclick=tog; ch.onkeydown=function(e){if(e.key=='Enter'||e.key==' '||e.key=='Spacebar'){e.preventDefault();tog();}};});
}

function renderTable(){
  var mets=R.metrics.filter(function(m){return !st.hidden[m.key];});
  function hsa(k){return ' tabindex="0" aria-sort="'+(st.sortKey==k?(st.asc?'ascending':'descending'):'none')+'"';}  // sortable-header a11y
  function sarr(k){return st.sortKey==k?(st.asc?icon('chevronUp','sort'):icon('chevronDown','sort')):'';}  // active-sort direction caret
  var head='<tr><th class="s" data-k="s"'+hsa('s')+'><input type="checkbox" id="cbAll" title="exclude all shown samples"><span class="hlab"> Sample</span>'+sarr('s')+'</th><th data-k="v"'+hsa('v')+'>QC'+sarr('v')+'</th>'+
    mets.map(function(m){var d=(R.defs[m.key]||[''])[0];return '<th data-k="'+m.key+'"'+hsa(m.key)+' title="'+esc(d)+'">'+esc(m.label)+(d?'<span class="infoi" title="'+esc(d)+'">i</span>':'')+sarr(m.key)+'</th>';}).join('')+
    '<th data-k="lineage"'+hsa('lineage')+' style="text-align:left">Lineage</th></tr>';
  // optional per-column filter row: one input per column (numeric ops on metric columns, substring otherwise)
  function cfIn(k,ph,lab,num){return '<input class="cfx" type="search" data-fk="'+k+'" value="'+esc(st.colf[k]||'')+'" placeholder="'+esc(ph)+'" aria-label="Filter '+esc(lab||k)+'"'+(num?' title="operators: &gt; &gt;= &lt; &lt;= = , a range 5-9 or 5..9; otherwise matches the text"':'')+'>';}
  var filtRow = st.showColF ? ('<tr class="colfilt"><th class="s">'+cfIn('s','name…','Sample')+'</th><th>'+cfIn('v','PASS/WARN…','QC status')+'</th>'+
    mets.map(function(m){return '<th>'+cfIn(m.key,'>50  5-9…',m.label,1)+'</th>';}).join('')+'<th>'+cfIn('lineage','L4…','Lineage')+'</th></tr>') : '';
  var linRank={}; (R.lineages||[]).forEach(function(l,i){linRank[l]=i;});
  function lr(s){return (s.lineage&&linRank[s.lineage]!=null)?linRank[s.lineage]:9999;}
  var shown=visible().filter(colMatch);
  var rows=shown.slice().sort(function(a,b){
    if(st.groupLin){var la=lr(a),lb=lr(b);if(la!=lb)return la-lb;}
    var k=st.sortKey,x=(k=='s')?a.s:(k=='v'?qcScore(a):a.m[k]),y=(k=='s')?b.s:(k=='v'?qcScore(b):b.m[k]),c;   // QC column sorts by severity, not the verdict string
    if(typeof x=='number'&&typeof y=='number')c=x-y;else c=String(x==null?'':x).localeCompare(String(y==null?'':y));return st.asc?c:-c;});
  var lastLin=null, ncol=mets.length+3;
  var TBL_CAP=400, total=rows.length, capped=total>TBL_CAP, draw=capped?rows.slice(0,TBL_CAP):rows;
  var body=draw.map(function(s){
    var pre='';
    if(st.groupLin){var lk=s.lineage||'NA'; if(lk!==lastLin){lastLin=lk;
      pre='<tr class="lingrp"><td class="s" colspan="'+ncol+'" style="text-align:left"><span class="ldot" style="background:'+linColor(s.lineage)+'"></span>'+esc(lk)+'</td></tr>';}}
    var badge=s.anc?'<span class="abadge" title="ancient (aDNA) sample">aDNA</span>':'';
    var tds='<td class="s" data-s="'+esc(s.s)+'"><input type="checkbox" class="cbx" data-s="'+esc(s.s)+'"'+(st.excl[s.s]?' checked':'')+' aria-label="basket '+esc(s.s)+'"><span class="sname" data-s="'+esc(s.s)+'" role="button" tabindex="0" aria-label="Open profile for '+esc(s.s)+'">'+esc(s.s)+'</span>'+badge+'</td><td data-v="'+s.v+'"><span class="v '+s.v+'">'+s.v+'</span></td>';
    mets.forEach(function(m){var v=s.m[m.key];
      if(v==null){tds+='<td class="na" data-v="">NA</td>';return;}
      var r=RANGES[m.key],nn=r[1]>r[0]?(v-r[0])/(r[1]-r[0]):0;nn=Math.max(0,Math.min(1,nn));var p=(nn*100).toFixed(1);
      tds+='<td data-v="'+v+'" style="background:linear-gradient(90deg,'+BAR[m.dir]+'2b 0 '+p+'%,#0000 '+p+'%)">'+fmt(v,m.kind)+'</td>';});
    tds+='<td data-v="'+esc(s.lineage||'')+'" style="text-align:left">'+(s.lineage?'<span class="ldot" style="background:'+linColor(s.lineage)+'"></span>':'')+esc(s.lineage||'NA')+'</td>';
    return pre+'<tr class="'+(st.hi==s.s?'hl':'')+'" data-s="'+esc(s.s)+'">'+tds+'</tr>';}).join('');
  var bodyOut=draw.length?body:'<tr><td colspan="'+ncol+'" style="text-align:left;color:#5f6f81;padding:14px 12px">No samples match the current filters.</td></tr>';
  var t=el('gstable'); t.innerHTML='<thead>'+head+filtRow+'</thead><tbody>'+bodyOut+'</tbody>';
  var af=anyFilterActive(), cntTxt=capped?('first '+TBL_CAP+' of '+total):(total+' / '+R.samples.length);
  el('nshown').innerHTML=cntTxt+' shown'+(af?' <a href="#" id="clrfilt" style="color:var(--accent);cursor:pointer;margin-left:7px;text-decoration:none">clear filters '+icon('x','sort')+'</a>':'');
  var cf=el('clrfilt'); if(cf)cf.onclick=function(e){e.preventDefault();clearAllFilters();};
  Array.prototype.forEach.call(t.querySelectorAll('th[data-k]'),function(th){function srt(){var k=th.getAttribute('data-k');if(st.sortKey==k)st.asc=!st.asc;else{st.sortKey=k;st.asc=(k=='s');}renderTable();saveState();}
    th.onclick=srt; th.onkeydown=function(e){if(e.key=='Enter'||e.key==' '||e.key=='Spacebar'){e.preventDefault();srt();}};});
  Array.prototype.forEach.call(t.querySelectorAll('.cfx'),function(inp){
    inp.onclick=function(e){e.stopPropagation();};
    inp.oninput=function(){var k=inp.getAttribute('data-fk'),pos=inp.selectionStart;st.colf[k]=inp.value;clearTimeout(_cfdb);_cfdb=setTimeout(function(){renderTable();
      var again=el('gstable').querySelector('.cfx[data-fk="'+k+'"]');if(again){again.focus();try{again.setSelectionRange(pos,pos);}catch(e){}}},140);};});
  Array.prototype.forEach.call(t.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(){setHi(tr.getAttribute('data-s'));};});
  // curation basket: exclusion checkboxes (stopPropagation so they don't sort/highlight)
  var cbAll=el('cbAll');
  if(cbAll){cbAll.checked=shown.length>0&&shown.every(function(s){return st.excl[s.s];});
    cbAll.onclick=function(e){e.stopPropagation();var on=cbAll.checked;shown.forEach(function(s){if(on)st.excl[s.s]=1;else delete st.excl[s.s];});renderTable();renderCuration();};}
  Array.prototype.forEach.call(t.querySelectorAll('.cbx'),function(cb){
    cb.onclick=function(e){e.stopPropagation();};
    cb.onchange=function(){var s=cb.getAttribute('data-s');if(cb.checked)st.excl[s]=1;else delete st.excl[s];renderCuration();
      var cba=el('cbAll');if(cba){cba.checked=shown.length>0&&shown.every(function(x){return st.excl[x.s];});}};});
  Array.prototype.forEach.call(t.querySelectorAll('.sname'),function(sp){sp.onclick=function(e){e.stopPropagation();openDetail(sp.getAttribute('data-s'));};
    sp.onkeydown=function(e){if(e.key=='Enter'||e.key==' '||e.key=='Spacebar'){e.preventDefault();e.stopPropagation();openDetail(sp.getAttribute('data-s'));}};});
}

var plotsZoom=30;   // beeswarm/bar/histogram plot height (px); the slider re-renders (dims are baked into the SVG)
function renderPlots(){
  var host=el('plots'); host.innerHTML='';
  var W=host.clientWidth||900, narrow=W<520, labelW=narrow?96:156, svgW=Math.max(narrow?190:240,W-labelW-4), padL=8, rightPad=64, pw=svgW-padL-rightPad, H=plotsZoom, cy=H/2;
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  DIST.forEach(function(pk){
    var mt=MET[pk]||{label:pk,kind:'float'};
    var rows=R.samples.filter(function(s){return s.m[pk]!=null;});
    var row=document.createElement('div'); row.className='bee';
    if(!rows.length){row.innerHTML='<span class="bl">'+esc(mt.label)+'</span><span class="nd">no data</span>';host.appendChild(row);return;}
    var lo=Math.min.apply(null,rows.map(function(s){return s.m[pk];})),hi=Math.max.apply(null,rows.map(function(s){return s.m[pk];}));
    if(hi<=lo)hi=lo+(Math.abs(lo)||1)*1e-3+1e-9;
    // acceptable-range band (moves with live thresholds); drawn behind marks. bar mode = ranked X, no band.
    var bandSVG='';(function(){var b=bandFor(pk,thr);if(!b||st.ptype=='bar')return;var blo=b[0],bhi=b[1];
      var x0=padL+(Math.max(blo,lo)-lo)/(hi-lo)*pw, x1=padL+(Math.min(bhi,hi)-lo)/(hi-lo)*pw, xw=x1-x0;
      if(xw>0.5){bandSVG='<rect x="'+x0.toFixed(1)+'" y="2" width="'+xw.toFixed(1)+'" height="'+(H-4)+'" fill="var(--bandfill)"/>';
        if(blo>lo)bandSVG+='<line x1="'+x0.toFixed(1)+'" y1="1" x2="'+x0.toFixed(1)+'" y2="'+(H-1)+'" stroke="var(--bandedge)" stroke-width="1"/>';
        if(bhi<hi)bandSVG+='<line x1="'+x1.toFixed(1)+'" y1="1" x2="'+x1.toFixed(1)+'" y2="'+(H-1)+'" stroke="var(--bandedge)" stroke-width="1"/>';}})();
    var inner='';
    if(st.ptype=='beeswarm'){
      inner=rows.map(function(s,i){var x=padL+(s.m[pk]-lo)/(hi-lo)*pw,j=((i*2654435761)%997)/997-0.5,y=cy+j*(H-9),
        big=(st.hi==s.s),dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],rr=big?4.7:(s.v!='PASS'?3.1:2.3),op=dim?0.1:(s.v!='PASS'?0.95:0.5),
        stk=(s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.3:0.6)+'"':'';
        return '<circle cx="'+x.toFixed(1)+'" cy="'+y.toFixed(1)+'" r="'+rr+'" fill="'+dotColor(s)+'" opacity="'+op+'"'+stk+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-pk="'+esc(pk)+'" data-val="'+s.m[pk]+'" data-lab="'+esc(mt.label)+'" data-kind="'+mt.kind+'"/>';}).join('');
      var mx=padL+(MED[pk]-lo)/(hi-lo)*pw;
      inner=bandSVG+'<line x1="'+padL+'" y1="'+cy+'" x2="'+(padL+pw)+'" y2="'+cy+'" stroke="'+TH.grid+'"/><line x1="'+mx.toFixed(1)+'" y1="4" x2="'+mx.toFixed(1)+'" y2="'+(H-4)+'" stroke="#64748b" stroke-dasharray="2 2"/>'+inner;
    }else if(st.ptype=='bar'){
      var sr=rows.slice().sort(function(a,b){return b.m[pk]-a.m[pk];}); var bw=pw/sr.length;
      inner=sr.map(function(s,i){var h=(s.m[pk]-Math.min(lo,0))/(hi-Math.min(lo,0))*(H-6),x=padL+i*bw,dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],big=(st.hi==s.s);
        return '<rect class="hit" x="'+x.toFixed(1)+'" y="'+(H-3-h).toFixed(1)+'" width="'+Math.max(bw-0.5,0.6).toFixed(1)+'" height="'+Math.max(h,0.5).toFixed(1)+'" fill="'+(big?''+TH.ink+'':dotColor(s))+'" opacity="'+(dim?0.12:(s.v!='PASS'?0.95:0.62))+'" data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-pk="'+esc(pk)+'" data-val="'+s.m[pk]+'" data-lab="'+esc(mt.label)+'" data-kind="'+mt.kind+'"/>';}).join('');
    }else{ // histogram
      var nb=Math.min(30,Math.max(8,Math.round(Math.sqrt(rows.length)))),cnt=zeros(nb);
      rows.forEach(function(s){var b=Math.floor((s.m[pk]-lo)/(hi-lo)*nb);if(b>=nb)b=nb-1;if(b<0)b=0;cnt[b]++;});
      var cm=Math.max.apply(null,cnt)||1,bw=pw/nb;
      inner=bandSVG+cnt.map(function(c,i){var h=c/cm*(H-6),x=padL+i*bw;return '<rect x="'+x.toFixed(1)+'" y="'+(H-3-h).toFixed(1)+'" width="'+Math.max(bw-1,0.6).toFixed(1)+'" height="'+Math.max(h,0.4).toFixed(1)+'" fill="'+BAR[mt.dir]+'" opacity="0.8"/>';}).join('');
    }
    var iq=IQR[pk],statTxt=iq?('med '+shortv(iq[1],mt.kind)+' · IQR '+shortv(iq[0],mt.kind)+' to '+shortv(iq[2],mt.kind)):'';
    row.innerHTML='<span class="bl">'+esc(mt.label)+(statTxt?'<span class="stat">'+esc(statTxt)+'</span>':'')+'</span><div class="plotwrap"><svg width="'+svgW+'" height="'+H+'">'+inner+
      '<text x="'+(padL+pw+6)+'" y="'+(cy+3)+'" font-size="9" fill="#94a3b8">'+shortv(hi,mt.kind)+'</text></svg></div>';
    host.appendChild(row);
  });
  host.style.setProperty('--beeh',(plotsZoom+10)+'px');
  var pz=el('plotszoom'); if(pz){ pz.value=plotsZoom; pz.oninput=function(){ plotsZoom=+this.value; renderPlots(); }; }
}

function renderScatter(){
  var host=el('scatter'); var cs=getComputedStyle(host); var avail=(host.clientWidth||560)-(parseFloat(cs.paddingLeft)||0)-(parseFloat(cs.paddingRight)||0); if(!(avail>0))avail=520; var cap=host.classList.contains('expanded')?880:700; var S=Math.max(240,Math.min(cap,avail)); var pad=42, plot=S-pad-14, H=(host.classList.contains('expanded')?Math.min(660,S):340), ph=H-pad-14;
  var xk=st.sx,yk=st.sy,xm=MET[xk],ym=MET[yk];
  var rows=R.samples.filter(function(s){return s.m[xk]!=null&&s.m[yk]!=null;});
  if(!rows.length){host.innerHTML='<div class="pad nd">no data for these axes</div>';return;}
  var xr=RANGES[xk],yr=RANGES[yk];
  function sx(v){return pad+(xr[1]>xr[0]?(v-xr[0])/(xr[1]-xr[0]):0.5)*plot;}
  function sy(v){return H-pad-(yr[1]>yr[0]?(v-yr[0])/(yr[1]-yr[0]):0.5)*ph;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  var dots=rows.map(function(s){var big=(st.hi==s.s),dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    return '<circle cx="'+sx(s.m[xk]).toFixed(1)+'" cy="'+sy(s.m[yk]).toFixed(1)+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="'+(dim?0.12:0.82)+'"'+((s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+s.m[xk]+'" data-y="'+s.m[yk]+'" data-xl="'+esc(xm.label)+'" data-yl="'+esc(ym.label)+'" data-xk="'+xm.kind+'" data-yk="'+ym.kind+'"/>';}).join('');
  var ticks='';[0,0.5,1].forEach(function(t){var gx=pad+t*plot,gy=H-pad-t*ph;
    ticks+='<line x1="'+gx+'" y1="'+pad+'" x2="'+gx+'" y2="'+(H-pad)+'" stroke="'+TH.grid+'"/><line x1="'+pad+'" y1="'+gy+'" x2="'+(pad+plot)+'" y2="'+gy+'" stroke="'+TH.grid+'"/>'+
    '<text x="'+gx+'" y="'+(H-pad+13)+'" font-size="9" fill="#94a3b8" text-anchor="middle">'+shortv(xr[0]+t*(xr[1]-xr[0]),xm.kind)+'</text>'+
    '<text x="'+(pad-6)+'" y="'+(gy+3)+'" font-size="9" fill="#94a3b8" text-anchor="end">'+shortv(yr[0]+t*(yr[1]-yr[0]),ym.kind)+'</text>';});
  // Spearman rho + two-sided p over the samples IN VIEW (matches the correlation-matrix scoping); a
  // quantitative read-out for any axis pair (e.g. a dose metric vs a QC metric).
  var cx=[],cy=[]; rows.forEach(function(s){if(vis[s.s]){cx.push(s.m[xk]);cy.push(s.m[yk]);}});
  var rho=spearman(cx,cy), pv=spearmanP(rho,cx.length);
  var corrCap=(cx.length>=4)?('Spearman &rho; = <b>'+(rho==null?'n/a':(rho>0?'':'−')+Math.abs(rho).toFixed(2))+'</b> &middot; p = '+pfmt(pv)+' &middot; n = '+cx.length+((pv!=null&&pv<0.05)?' <span class="sc-sig">significant</span>':''))
    :('n = '+cx.length+' — need &ge; 4 samples in view for a correlation');
  host.innerHTML='<svg width="'+S+'" height="'+H+'" id="scsvg" style="display:block;max-width:100%;margin:0 auto">'+
    '<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(pad+plot)+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/><line x1="'+pad+'" y1="'+pad+'" x2="'+pad+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/>'+
    ticks+dots+
    '<text x="'+(pad+plot/2)+'" y="'+(H-6)+'" font-size="11" fill="'+TH.mut+'" text-anchor="middle">'+esc(xm.label)+'</text>'+
    '<text x="12" y="'+(pad+ph/2)+'" font-size="11" fill="'+TH.mut+'" text-anchor="middle" transform="rotate(-90 12 '+(pad+ph/2)+')">'+esc(ym.label)+'</text></svg>'+
    '<div class="sc-corr" title="Spearman rank correlation between the two selected axes over the samples currently in view, with a two-sided p-value (Student-t approximation) and the sample count. Descriptive; not corrected for multiple comparisons.">'+corrCap+'</div>'+
    colorLegend();
  SGEO={pad:pad,plot:plot,ph:ph,H:H,xr:xr,yr:yr,xk:xk,yk:yk};   // for the rubber-band select inverse-mapping
  var scsvg=el('scsvg'); if(scsvg){var ov=document.createElementNS('http://www.w3.org/2000/svg','rect');
    ov.setAttribute('id','scbrush');ov.setAttribute('fill','rgba(14,139,168,.12)');ov.setAttribute('stroke','#0e8ba8');
    ov.setAttribute('stroke-dasharray','3 2');ov.setAttribute('pointer-events','none');ov.style.display='none';scsvg.appendChild(ov);}
}

// ---- Spearman correlation matrix across the core metrics (click a cell -> loads that pair into the scatter) ----
function rankvec(a){var idx=a.map(function(v,i){return [v,i];}).sort(function(x,y){return x[0]-y[0];});
  var r=new Array(a.length),k=0;
  while(k<idx.length){var j=k;while(j+1<idx.length&&idx[j+1][0]===idx[k][0])j++;
    var avg=(k+j)/2+1;for(var t=k;t<=j;t++)r[idx[t][1]]=avg;k=j+1;}
  return r;}
function spearman(x,y){var n=x.length; if(n<4)return null;
  var rx=rankvec(x),ry=rankvec(y),mx=0,my=0,i;
  for(i=0;i<n;i++){mx+=rx[i];my+=ry[i];} mx/=n;my/=n;
  var sxy=0,sxx=0,syy=0;
  for(i=0;i<n;i++){var dx=rx[i]-mx,dy=ry[i]-my;sxy+=dx*dy;sxx+=dx*dx;syy+=dy*dy;}
  return (sxx>0&&syy>0)?sxy/Math.sqrt(sxx*syy):null;}
function corrCol(r){if(r==null)return TH.cellnull;var a=Math.abs(r),base=r>=0?[224,84,79]:[79,131,194],w=isDark()?[30,42,56]:[247,249,252];
  return 'rgb('+w.map(function(c,i){return Math.round(c+(base[i]-c)*a);}).join(',')+')';}
function renderCorr(){
  var host=el('corr_body'); if(!host)return;
  var keys=DIST.filter(function(k){return MET[k];}).slice(0,16);
  var pool=visible();
  if(pool.length<4||keys.length<2){host.innerHTML='<div class="pad nd">not enough data for a correlation matrix (need >= 4 samples in view).</div>';return;}
  var vals={}; keys.forEach(function(k){vals[k]=pool.map(function(s){return s.m[k];});});
  var n=keys.length, cell=Math.max(16,Math.min(34,Math.floor((Math.min(host.clientWidth||560,host.classList.contains('expanded')?1100:940)-110)/n)));
  var padL=96,padT=8, W=padL+n*cell+8, H=padT+n*cell+128;
  var svg='<svg width="'+W+'" height="'+H+'" id="corrsvg" style="display:block">';   // no max-width: let #corr_body{overflow-x:auto} scroll instead of clipping the right columns
  keys.forEach(function(k,j){var cx=padL+j*cell+cell/2;
    svg+='<text x="'+cx+'" y="'+(padT+n*cell+12)+'" font-size="8.5" fill="'+TH.mut+'" text-anchor="end" transform="rotate(-55 '+cx+' '+(padT+n*cell+12)+')">'+esc(MET[k].label)+'</text>';});
  keys.forEach(function(k,i){var cy=padT+i*cell+cell/2;
    svg+='<text x="'+(padL-6)+'" y="'+(cy+3)+'" font-size="8.5" fill="'+TH.mut+'" text-anchor="end">'+esc(MET[k].label)+'</text>';});
  for(var i=0;i<n;i++)for(var j=0;j<n;j++){
    var cx=padL+j*cell,cy=padT+i*cell,r;
    if(i===j){r=vals[keys[i]].some(function(v){return v!=null;})?1:null;}   // grey out an all-NA metric
    else{var xs=[],ys=[]; for(var t=0;t<pool.length;t++){var xv=vals[keys[j]][t],yv=vals[keys[i]][t];
        if(xv!=null&&yv!=null){xs.push(xv);ys.push(yv);}} r=spearman(xs,ys);}
    svg+='<rect x="'+cx+'" y="'+cy+'" width="'+(cell-1)+'" height="'+(cell-1)+'" rx="2" fill="'+corrCol(r)+'"'+
      ' data-xk="'+keys[j]+'" data-yk="'+keys[i]+'" data-r="'+(r==null?'':r.toFixed(2))+'"'+(i!==j?' style="cursor:pointer"':'')+'/>';
    if(cell>=22&&r!=null)svg+='<text x="'+(cx+(cell-1)/2)+'" y="'+(cy+(cell-1)/2+3)+'" font-size="7.5" fill="'+(Math.abs(r)>0.55?'#fff':'#5b6b7e')+'" text-anchor="middle" pointer-events="none">'+(r>0?'':'-')+Math.abs(r).toFixed(1).replace('0.','.')+'</text>';
  }
  var ly=padT+n*cell+94;
  svg+='<text x="'+padL+'" y="'+(ly-4)+'" font-size="8.5" fill="#94a3b8">Spearman rho</text>';
  for(var g=0;g<=20;g++){var rr=-1+g/10; svg+='<rect x="'+(padL+g*7)+'" y="'+ly+'" width="7" height="9" fill="'+corrCol(rr)+'"/>';}
  svg+='<text x="'+padL+'" y="'+(ly+20)+'" font-size="8" fill="#94a3b8">-1</text>'+
       '<text x="'+(padL+70)+'" y="'+(ly+20)+'" font-size="8" fill="#94a3b8" text-anchor="middle">0</text>'+
       '<text x="'+(padL+140)+'" y="'+(ly+20)+'" font-size="8" fill="#94a3b8" text-anchor="end">+1</text>';
  host.innerHTML=svg+'</svg>';
}

var consZoom=26;   // consensus per-sample bar row height (px); CSS var --stackh, no re-render
function renderStacks(){
  var host=el('stacks');
  var rows=R.samples.filter(function(s){return s.m.callable_pct!=null||s.m.missing_pct!=null;})
    .slice().sort(function(a,b){return (b.m.missing_pct||0)-(a.m.missing_pct||0);});
  if(!rows.length){host.innerHTML='<div class="pad nd">no consensus data</div>';return;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  host.innerHTML=rows.map(function(s){var cal=s.m.callable_pct||0,iup=s.m.iupac_pct||0,mis=Math.max(0,100-cal-iup);
    var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    return '<div class="stack" data-s="'+esc(s.s)+'" style="opacity:'+(dim?0.25:1)+(st.hi==s.s?';background:'+TH.hl:'')+'">'+
     '<span class="sl">'+esc(s.s)+'</span><div class="sb">'+
     '<div style="width:'+cal.toFixed(2)+'%;background:#22a06b" title="callable"></div>'+
     '<div style="width:'+iup.toFixed(2)+'%;background:#e6b25a" title="IUPAC"></div>'+
     '<div style="width:'+mis.toFixed(2)+'%;background:#cbd5e1" title="missing"></div></div>'+
     '<span class="sv">'+mis.toFixed(1)+'%</span></div>';}).join('');
  Array.prototype.forEach.call(host.querySelectorAll('.stack'),function(d){d.onclick=function(){setHi(d.getAttribute('data-s'));};});
  function applyZ(){ host.style.setProperty('--stackh',consZoom+'px'); host.style.setProperty('--sbh',Math.round(consZoom*0.5)+'px'); }
  applyZ();
  var cz=el('conszoom'); if(cz){ cz.value=consZoom; cz.oninput=function(){ consZoom=+this.value; applyZ(); }; }
}

function qcScore(s){var fails=s.f.filter(function(f){return FAILF[f];}).length;return fails*100+(s.f.length-fails)*10;}
function renderFlags(){
  var fl=R.samples.filter(function(s){return s.v!='PASS';}).sort(function(a,b){return qcScore(b)-qcScore(a)||a.s.localeCompare(b.s);});
  var body=fl.length?fl.map(function(s){return '<tr data-s="'+esc(s.s)+'">'+
    '<td class="s">'+esc(s.s)+(s.anc?'<span class="abadge">aDNA</span>':'')+'</td><td><span class="v '+s.v+'">'+s.v+'</span></td>'+
    '<td class="flags"><div class="fchips">'+s.f.map(function(f){return '<span class="chip'+(FAILF[f]?' failc':'')+'" data-f="'+f+'" title="'+esc(flagWhy(s,f))+'" style="cursor:pointer">'+f+'</span>';}).join(' ')+'</div>'+
      '<div class="flagrsn">'+s.f.map(function(f){return '<span>'+esc(flagWhy(s,f))+'</span>';}).join('')+'</div></td></tr>';}).join('')
    :'<tr><td colspan="3" style="text-align:left;color:#16a34a;padding:10px">All samples pass at the current thresholds.</td></tr>';
  var t=el('flagtable');
  t.innerHTML='<thead><tr><th class="s">Sample (worst first)</th><th>QC</th><th style="text-align:left">Flags &amp; reason</th></tr></thead><tbody>'+body+'</tbody>';
  Array.prototype.forEach.call(t.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(e){
    if(e.target.classList.contains('chip')){var f=e.target.getAttribute('data-f');st.flagFilter=(st.flagFilter==f?null:f);renderAll();el('gstats').scrollIntoView();return;}
    setHi(tr.getAttribute('data-s'));el('gstats').scrollIntoView();};});
  el('nflag').textContent=fl.length;
  var bw=el('basketFlagged'); if(bw){bw.onclick=function(){fl.forEach(function(s){st.excl[s.s]=1;});renderTable();renderCuration();};}
}

// ---- curation basket + exclusion exports ----
function nExcl(){return R.samples.filter(function(s){return st.excl[s.s];}).length;}
