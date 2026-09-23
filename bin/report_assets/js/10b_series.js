// ---- Series: what each group gained since its first time point ----
// A group of the samplesheet (a patient, a passage line) followed over time. Each later sample's
// fixed SNPs are read against the group's first time point, where the SNP matrix says whether a
// site was read without the allele (new) or not read at all (unknown).
var SERPAL=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#e69f00','#0072b2','#d55e00','#009e73','#cc79a7','#882255'];
function serColorBy(){var S=R.series,m={},k=0;if(!S)return function(){return '#4f83c2';};
  var key=function(g){return S.tx_field?(g.tx||'NA'):g.group;};
  S.groups.map(key).sort().forEach(function(v){if(!(v in m)){m[v]=SERPAL[k%SERPAL.length];k++;}});
  var f=function(g){return m[key(g)];};f.map=m;return f;}
function serGained(r){return r['new'].length+r.risen.length;}
// A sample the QC fails or places outside its series (another lineage, far from its group) is not
// the series evolving: it is drawn hollow and left out of the line and of the scale, so one
// mislabelled sample carrying 2,000 SNPs does not flatten every real series to the floor.
var SER_OUT={GROUP_MISMATCH:1,LINEAGE_MISMATCH:1};
function serOutside(sid){var s=null;R.samples.forEach(function(x){if(x.s===sid)s=x;});
  return !!s&&(s.v==='FAIL'||(s.f||[]).some(function(f){return SER_OUT[f];}));}
function serMutLabel(site){var x=(R.series.sites||{})[site]||{};return (x.gene?x.gene:site)+(x.aa?' '+x.aa.replace(/^p\./,''):'');}
function renderSeries(){
  var host=el('gains_body'),S=R.series; if(!host)return; if(!(S&&S.groups&&S.groups.length)){host.innerHTML='';return;}
  var col=serColorBy(),numeric=S.groups.every(function(g){return g.tnum0!=null&&g.rows.every(function(r){return r.tnum!=null;});});
  var pts=[];S.groups.forEach(function(g){pts.push(numeric?g.tnum0:0);g.rows.forEach(function(r,i){pts.push(numeric?r.tnum:i+1);});});
  var x0=Math.min.apply(null,pts),x1=Math.max.apply(null,pts),ymax=1,nOut=0;
  S.groups.forEach(function(g){g.rows.forEach(function(r){if(serOutside(r.s))nOut++;else ymax=Math.max(ymax,serGained(r));});});
  var W=Math.max(320,Math.min((host.clientWidth||900)-24,980)),H=260,pl=40,pr=12,pt=12,pb=28,iw=W-pl-pr,ih=H-pt-pb;
  function X(v){return pl+(x1>x0?(v-x0)/(x1-x0):0.5)*iw;} function Y(v){return pt+ih-v/ymax*ih;}
  var svg='<svg width="'+W+'" height="'+H+'">';
  [0,0.5,1].forEach(function(t){var v=Math.round(t*ymax);svg+='<line x1="'+pl+'" x2="'+(W-pr)+'" y1="'+Y(v).toFixed(1)+'" y2="'+Y(v).toFixed(1)+'" stroke="'+TH.grid+'"/>'+
    '<text x="'+(pl-6)+'" y="'+(Y(v)+3).toFixed(1)+'" font-size="9" fill="'+TH.mut+'" text-anchor="end">'+v+'</text>';});
  var ticks={};pts.forEach(function(v){ticks[v]=1;});
  Object.keys(ticks).map(Number).sort(function(a,b){return a-b;}).forEach(function(v){
    svg+='<text x="'+X(v).toFixed(1)+'" y="'+(H-10)+'" font-size="9" fill="'+TH.mut+'" text-anchor="middle">'+(numeric?v:'')+'</text>';});
  S.groups.forEach(function(g){var c=col(g),byT={},order=[];
    // the line runs through the median of the samples counted at each time point
    g.rows.forEach(function(r,i){var x=numeric?r.tnum:i+1;if(serOutside(r.s))return;if(!(x in byT)){byT[x]=[];order.push(x);}byT[x].push(serGained(r));});
    var line=[[numeric?g.tnum0:0,0]].concat(order.sort(function(a,b){return a-b;}).map(function(x){var v=byT[x].slice().sort(function(a,b){return a-b;});return [x,v[Math.floor(v.length/2)]];}));
    svg+='<polyline points="'+line.map(function(p){return X(p[0]).toFixed(1)+','+Y(p[1]).toFixed(1);}).join(' ')+'" fill="none" stroke="'+c+'" stroke-width="1.6" opacity="0.85"/>';
    svg+='<circle cx="'+X(line[0][0]).toFixed(1)+'" cy="'+Y(0).toFixed(1)+'" r="3.2" fill="'+c+'" data-tip="'+(esc(g.group)+' &#183; start ('+esc(g.t0)+')').replace(/"/g,'&quot;')+'"/>';
    g.rows.forEach(function(r,i){var x=numeric?r.tnum:i+1,out=serOutside(r.s),v=serGained(r),y=out?Math.min(v,ymax):v;
      var lab=esc(r.s)+' &#183; '+esc(r.time)+': '+r['new'].length+' new, '+r.risen.length+' risen'+(r.unknown?', '+r.unknown+' unknown':'')+(r.lost.length?', '+r.lost.length+' lost':'')+(out?' &#183; not counted: the QC fails it or places it outside its series':'');
      svg+='<circle cx="'+X(x).toFixed(1)+'" cy="'+Y(y).toFixed(1)+'" r="3.2" '+(out?'fill="none" stroke="'+c+'" stroke-width="1.4"':'fill="'+c+'"')+' data-tip="'+lab.replace(/"/g,'&quot;')+'"/>';});});
  svg+='<text x="'+pl+'" y="'+(pt+2)+'" font-size="9" fill="'+TH.mut+'">fixed SNPs gained since the first time point</text></svg>';
  var legend='<div class="legend">'+Object.keys(col.map).map(function(k){return '<span><i style="background:'+col.map[k]+'"></i>'+esc(k)+'</span>';}).join('')+
    (nOut?'<span><i style="background:transparent;border:1.4px solid '+TH.mut+'"></i>'+nOut+' not counted (FAIL, or outside its series), capped at the top of the scale</span>':'')+
    '<span style="margin-left:auto">'+(S.tx_field?'coloured by '+esc(S.tx_field):'coloured by group')+'; the line is the median per time point; hover a point for its sample</span></div>';
  var rows=[];S.groups.forEach(function(g){g.rows.forEach(function(r,i){
    var muts=r['new'].map(function(s){return '<b>'+esc(serMutLabel(s))+'</b>';}).concat(r.risen.map(function(s){return esc(serMutLabel(s))+' <span class="na">(risen)</span>';}));
    rows.push('<tr><td class="s" style="text-align:left">'+(i===0?'<i class="ser-dot" style="background:'+col(g)+'"></i>'+esc(g.group):'')+'</td>'+
      (S.tx_field?'<td style="text-align:left">'+(i===0?esc(g.tx||''):'')+'</td>':'')+
      '<td style="text-align:left">'+esc(r.s)+(serOutside(r.s)?' <span class="na" title="the QC fails it or places it outside its series: not counted in the chart">(not counted)</span>':'')+'</td><td>'+esc(r.time)+'</td><td>'+r['new'].length+'</td><td>'+r.risen.length+'</td>'+
      '<td'+(r.unknown?' title="the first time point was not read at these sites"':'')+'>'+r.unknown+'</td>'+
      '<td'+(r.lost.length?' style="color:var(--warn);font-weight:600" title="'+esc(r.lost.map(serMutLabel).join(', '))+'"':'')+'>'+r.lost.length+'</td>'+
      '<td style="text-align:left">'+(muts.length?_list(muts,8):'<span class="na">none</span>')+'</td></tr>');});});
  host.innerHTML='<div class="hot-note" style="margin:0 0 8px">Fixed SNPs (allele fraction '+S.fixed+' or more) each sample carries that the first time point of its series did not: '+
    '<b>new</b> where that first time point was read at the site without the allele, <b>risen</b> where it held the allele as a minority, <b>unknown</b> where it was not read there. '+
    '<b>Lost</b> counts the SNPs fixed at the start and read without the allele later: a reversion, or a sample that does not belong to its series. '+
    (S.checked?'':'<b>The SNP matrix was not available</b>, so a site not called at the first time point could not be told from one not read there, and counts as new. ')+
    'SNPs only; indels and deletions are not counted.</div>'+
    '<div class="ser-chart">'+svg+'</div>'+legend+
    '<div class="gtable" style="max-height:46vh;margin-top:10px"><table><thead><tr><th class="s" style="text-align:left">Series</th>'+(S.tx_field?'<th style="text-align:left">'+esc(S.tx_field)+'</th>':'')+
    '<th style="text-align:left">Sample</th><th>Time</th><th>New</th><th>Risen</th><th>Unknown</th><th>Lost</th><th style="text-align:left">Gained (new in bold)</th></tr></thead><tbody>'+rows.join('')+'</tbody></table></div>';
  Array.prototype.forEach.call(host.querySelectorAll('circle[data-tip]'),function(c){
    c.onmousemove=function(e){tip(c.getAttribute('data-tip'),e.clientX,e.clientY);};c.onmouseleave=function(){tip('');};});
}
// Resistance along each series: the WHO grade 1-2 mutations at every time point, those acquired
// since the first one set apart. Lineage markers are there from the start, so they read as inherited.
function renderDrSeries(){
  var host=el('drseries_body'),S=R.series,D=R.dr; if(!host)return;
  if(!(S&&S.groups&&D&&D.calls)){host.innerHTML='';return;}
  var by={};D.calls.forEach(function(c){if(c.gn==null||c.gn>2)return;(by[c.s]=by[c.s]||{})[c.gene+' '+c.mutation]=c.drug;});
  var has={};D.calls.forEach(function(c){has[c.s]=1;});
  var out=[];
  S.groups.forEach(function(g){
    var base={};g.base.forEach(function(b){var m=by[b]||{};Object.keys(m).forEach(function(k){base[k]=m[k];});});
    var any=Object.keys(base).length>0,cells=[];
    g.rows.forEach(function(r){var m=by[r.s]||{},keys=Object.keys(m),acq=keys.filter(function(k){return !(k in base);});
      var lost=has[r.s]?Object.keys(base).filter(function(k){return !(k in m);}):[];
      if(keys.length||lost.length)any=true;
      cells.push({r:r,acq:acq,kept:keys.filter(function(k){return k in base;}),lost:lost,m:m});});
    if(any)out.push({g:g,base:base,cells:cells});});
  if(!out.length){host.innerHTML='<div class="hot-note">No series carries a grade 1&#8211;2 resistance mutation at any time point.</div>';return;}
  function mut(k,drug,cls){return '<span class="ser-mut'+(cls?' '+cls:'')+'" title="'+esc(drug||'')+'">'+esc(k)+(drug?' <span class="na">'+esc(drug)+'</span>':'')+'</span>';}
  var nAcq=0;out.forEach(function(o){o.cells.forEach(function(c){nAcq+=c.acq.length?1:0;});});
  host.innerHTML='<div class="hot-note" style="margin:0 0 8px">Grade 1&#8211;2 mutations of the WHO catalogue at each time point of each series. <b>Acquired</b> ones were not there at the first time point; a crossed-out one was there and is no longer called. '+
    nAcq+' later '+_plural(nAcq,'sample carries','samples carry')+' at least one acquired mutation. A genomic screen, not a susceptibility result.</div>'+
    '<div class="gtable" style="max-height:60vh"><table class="ser-tbl"><thead><tr><th class="s" style="text-align:left">Series</th>'+(S.tx_field?'<th style="text-align:left">'+esc(S.tx_field)+'</th>':'')+'<th style="text-align:left">First time point</th><th style="text-align:left">Later time points</th></tr></thead><tbody>'+
    out.map(function(o){var bk=Object.keys(o.base);
      return '<tr><td class="s" style="text-align:left">'+esc(o.g.group)+'</td>'+(S.tx_field?'<td style="text-align:left">'+esc(o.g.tx||'')+'</td>':'')+
        '<td style="text-align:left"><div class="na">'+esc(o.g.t0)+'</div>'+(bk.length?bk.map(function(k){return mut(k,o.base[k],'');}).join(' '):'<span class="na">none</span>')+'</td>'+
        '<td style="text-align:left">'+o.cells.map(function(c){return '<div class="ser-tp"><span class="na">'+esc(c.r.time)+' &#183; '+esc(c.r.s)+'</span> '+
          c.acq.map(function(k){return mut(k,c.m[k],'acq');}).concat(c.kept.map(function(k){return mut(k,c.m[k],'');})).concat(c.lost.map(function(k){return mut(k,o.base[k],'lost');})).join(' ')+
          (c.acq.length+c.kept.length+c.lost.length?'':'<span class="na">'+(has[c.r.s]?'none':'no resistance calls for this sample')+'</span>')+'</div>';}).join('')+'</td></tr>';}).join('')+
    '</tbody></table></div>';
}
