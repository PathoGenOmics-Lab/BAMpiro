// ---- Relatedness: SNP distances between the consensus sequences, clustered at a threshold ----
// The distances count only positions both samples called, so a gap is never a difference. A pair
// that compared too few of its reference's variable positions says little either way: it neither
// links a cluster nor orders the heatmap.
var REL={ref:null,thr:null};
var RELPAL=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#e69f00','#0072b2','#d55e00','#009e73','#cc79a7','#882255'];
function relSec(){return R.relatedness||null;}
function relRef(){var S=relSec();if(!S)return null;var ks=Object.keys(S.refs||{});if(!ks.length)return null;
  if(!REL.ref||!S.refs[REL.ref])REL.ref=ks.slice().sort(function(a,b){return S.refs[b].samples.length-S.refs[a].samples.length;})[0];
  return S.refs[REL.ref];}
function relThr(){var S=relSec();if(REL.thr==null)REL.thr=(S&&S.threshold!=null)?S.threshold:12;return REL.thr;}
function relMinCmp(){var S=relSec();return 100*((S&&S.min_compared!=null)?S.min_compared:0.5);}
function relIdx(n,i,j){if(i>j){var t=i;i=j;j=t;}return i*n-i*(i+1)/2+(j-i-1);}   // upper triangle, row by row
function relD(r,i,j){return i==j?0:r.snps[relIdx(r.samples.length,i,j)];}
function relC(r,i,j){return i==j?100:r.cmp[relIdx(r.samples.length,i,j)];}
function relOk(r,i,j,minCmp){var c=relC(r,i,j);return c!=null&&c>=minCmp;}
// Single-linkage clusters at the threshold (union-find), largest first.
function relClusters(r,thr,minCmp){var n=r.samples.length,p=[],i,j;for(i=0;i<n;i++)p.push(i);
  function f(x){while(p[x]!==x){p[x]=p[p[x]];x=p[x];}return x;}
  for(i=0;i<n;i++)for(j=i+1;j<n;j++){var d=relD(r,i,j);if(d!=null&&d<=thr&&relOk(r,i,j,minCmp)){var a=f(i),b=f(j);if(a!==b)p[a]=b;}}
  var by={},keys=[];for(i=0;i<n;i++){var k=f(i);if(!by[k]){by[k]=[];keys.push(k);}by[k].push(i);}
  return keys.map(function(k){return by[k];}).sort(function(a,b){return (b.length-a.length)||(a[0]-b[0]);});}
// A leaf order in which every single-linkage cluster is contiguous at any threshold: merge the
// closest pairs first (Kruskal) and join the two member lists at each merge.
function relOrder(r,minCmp){var n=r.samples.length,e=[],i,j;
  for(i=0;i<n;i++)for(j=i+1;j<n;j++){var d=relD(r,i,j);if(d!=null&&relOk(r,i,j,minCmp))e.push([d,i,j]);}
  e.sort(function(a,b){return a[0]-b[0];});
  var p=[],mem=[];for(i=0;i<n;i++){p.push(i);mem.push([i]);}
  function f(x){while(p[x]!==x){p[x]=p[p[x]];x=p[x];}return x;}
  e.forEach(function(x){var a=f(x[1]),b=f(x[2]);if(a===b)return;p[a]=b;mem[b]=mem[b].concat(mem[a]);mem[a]=null;});
  var out=[];for(i=0;i<n;i++)if(p[i]===i&&mem[i])out=out.concat(mem[i]);
  return out;}
function relGroupColor(){var S=relSec(),g=(S&&S.group)||{},m={},k=0;
  Object.keys(g).map(function(s){return g[s];}).sort().forEach(function(v){if(!(v in m)){m[v]=RELPAL[k%RELPAL.length];k++;}});
  return function(s){var v=g[s];return v?m[v]:null;};}
// Within the threshold: teal, darkest where identical. Beyond it: grey, darker as the distance grows (log).
function relColor(d,thr,max){if(d==null)return TH.cellnull;
  if(d<=thr){var t=thr>0?d/thr:0;return 'rgb('+Math.round(12+t*110)+','+Math.round(110+t*85)+','+Math.round(132+t*70)+')';}
  var lt=Math.log(d-thr+1)/Math.log(Math.max(2,max-thr+1)),g=isDark()?Math.round(60+lt*110):Math.round(232-lt*140);
  return 'rgb('+g+','+g+','+Math.min(255,g+8)+')';}
function relHistogram(r,thr,minCmp){
  // The threshold is one of the bin edges, so every bar is wholly within it or wholly beyond.
  var edges=[0,1,2,3,6,13,26,51,101,501,1001].filter(function(e){return e!==thr+1;}).concat([thr+1]).sort(function(a,b){return a-b;});
  edges.push(Infinity);
  var labels=edges.slice(0,-1).map(function(e,k){var hi=edges[k+1];
    return hi===Infinity?('> '+(e-1).toLocaleString('en-US')):(hi-1===e?String(e):e.toLocaleString('en-US')+'-'+(hi-1).toLocaleString('en-US'));});
  var cnt=zeros(labels.length),n=r.samples.length,tot=0,i,j,k;
  for(i=0;i<n;i++)for(j=i+1;j<n;j++){var d=relD(r,i,j);if(d==null||!relOk(r,i,j,minCmp))continue;tot++;
    for(k=0;k<labels.length;k++)if(d>=edges[k]&&d<edges[k+1]){cnt[k]++;break;}}
  var mx=Math.max.apply(null,cnt.concat([1]));
  return '<div class="rel-hist">'+labels.map(function(l,k){var within=edges[k+1]-1<=thr;
    return '<div class="rel-bar" title="'+cnt[k]+' pair(s) at '+l+' SNPs"><span class="rel-bv">'+(cnt[k]||'')+'</span>'+
      '<i style="height:'+Math.round(4+cnt[k]/mx*70)+'px;background:'+(within?'#0e8ba8':(isDark()?'#4b5b6e':'#b8c2cf'))+'"></i><span class="rel-bl">'+l+'</span></div>';}).join('')+
    '</div><div class="hot-note">'+tot.toLocaleString('en-US')+' pairs, in SNPs; teal = within the threshold of '+thr+'.</div>';}
// The controls are built once: rebuilding them on every redraw replaced the slider under the
// pointer mid-drag and dropped the keyboard focus. The slider's range is fixed for the same reason.
function relControls(host,S){
  if(el('relthr'))return;
  var refs=Object.keys(S.refs),smax=Math.max(50,Math.min(500,4*(S.threshold||12)));
  host.innerHTML='<div class="rel-ctl">'+(refs.length>1?('<span class="seg" id="relrefs">'+refs.map(function(k){return '<button data-r="'+esc(k)+'"'+(k==REL.ref?' class="on"':'')+'>'+esc(k)+' ('+S.refs[k].samples.length+')</button>';}).join('')+'</span>'):'')+
    '<label class="rel-thr" title="SNPs within which two samples are drawn as one cluster">cluster within <input type="number" id="relthr" min="0" max="100000" step="1" value="'+relThr()+'"> SNPs'+
    '<input type="range" id="relthrr" min="0" max="'+smax+'" step="1" value="'+Math.min(relThr(),smax)+'"></label></div><div id="rel_content"></div>';
  function setThr(v,from){v=Math.max(0,Math.round(+v||0));REL.thr=v;
    if(from!=='num')el('relthr').value=v; if(from!=='range')el('relthrr').value=Math.min(v,smax);
    renderRelContent();renderRelClusters();}
  el('relthr').onchange=function(){setThr(this.value,'num');};
  el('relthrr').oninput=function(){var v=this.value;clearTimeout(REL._db);REL._db=setTimeout(function(){setThr(v,'range');},120);};
  Array.prototype.forEach.call(host.querySelectorAll('#relrefs button'),function(b){b.onclick=function(){REL.ref=b.getAttribute('data-r');
    Array.prototype.forEach.call(host.querySelectorAll('#relrefs button'),function(x){x.classList.toggle('on',x==b);});
    renderRelContent();renderRelClusters();};});
}
function renderRelatedness(){
  var host=el('rel_body'),S=relSec(); if(!host)return;
  if(!(S&&relRef())){host.innerHTML='';return;}
  relControls(host,S); renderRelContent(); renderRelClusters();
}
function renderRelContent(){
  var box=el('rel_content'),S=relSec(),r=relRef(); if(!(box&&S&&r))return;
  // Drawn only once it can be seen: the load-time pass over every page would otherwise paint a
  // thousand-sample heatmap into a hidden page, and opening the page would paint it again.
  if(box.offsetParent===null&&!document.body.classList.contains('print-all'))return;
  var thr=relThr(),minCmp=relMinCmp(),n=r.samples.length,key=REL.ref+'|'+minCmp;
  // The order and the largest distance do not depend on the threshold: worked out once per reference.
  REL.cache=REL.cache||{};
  var C=REL.cache[key]||(REL.cache[key]={order:relOrder(r,minCmp),max:(function(){var m=0,i,j;for(i=0;i<n;i++)for(j=i+1;j<n;j++){var d=relD(r,i,j);if(d!=null&&d>m)m=d;}return m;})()});
  var order=C.order,max=C.max,gcol=relGroupColor();
  var cl=relClusters(r,thr,minCmp),multi=cl.filter(function(c){return c.length>1;}),alone=cl.filter(function(c){return c.length==1;}).length;
  box.innerHTML='<div class="hot-note" style="margin:4px 0 8px"><b>'+multi.length+'</b> '+_plural(multi.length,'cluster')+' of two or more samples within '+thr+' SNPs'+
    (multi.length?(', the largest of '+multi[0].length):'')+'; '+alone+' '+_plural(alone,'sample stands','samples stand')+' alone. '+n+' samples on '+esc(REL.ref)+
    ', '+r.variable.toLocaleString('en-US')+' variable positions. Pairs that both called under '+minCmp+'% of them are left out of the clusters and shown blank.</div>'+
    '<div class="rel-wrap"><canvas id="relcanvas"></canvas></div>'+
    '<div class="legend"><span><i style="background:'+relColor(0,thr,max)+'"></i>identical</span><span><i style="background:'+relColor(thr,thr,max)+'"></i>'+thr+' SNPs</span>'+
    '<span><i style="background:'+relColor(max,thr,max)+'"></i>'+max.toLocaleString('en-US')+' SNPs</span><span><i style="background:'+TH.cellnull+'"></i>too few positions compared</span>'+
    (Object.keys(S.group||{}).length?'<span>left strip = group</span>':'')+'<span style="margin-left:auto">samples in single-linkage order; hover a cell for the pair</span></div>'+
    relHistogram(r,thr,minCmp);
  // A square that fits the screen both ways, so the last rows are never scrolled out of sight.
  var cv=el('relcanvas'),W=Math.max(240,Math.min((box.clientWidth||900)-24,980)),H=Math.max(240,Math.min(760,Math.round((window.innerHeight||900)*0.72))),
      strip=10,cs=Math.max(1,Math.floor(Math.min(W-strip-4,H)/Math.max(1,n))),side=cs*n;
  var dpr=window.devicePixelRatio||1; cv.width=(side+strip+4)*dpr; cv.height=side*dpr; cv.style.width=(side+strip+4)+'px'; cv.style.height=side+'px';
  var ctx=cv.getContext&&cv.getContext('2d');
  if(ctx){ctx.scale(dpr,dpr);
    var cc={};   // one colour string per distinct distance, not one per cell
    for(var a=0;a<n;a++){var gc=gcol(r.samples[order[a]]);ctx.fillStyle=gc||TH.cellnull;ctx.fillRect(0,a*cs,strip-2,cs);
      for(var b=0;b<n;b++){var ia=order[a],ib=order[b],d=relD(r,ia,ib),okp=ia==ib||relOk(r,ia,ib,minCmp);
        ctx.fillStyle=okp?(d==null?TH.cellnull:(cc[d]||(cc[d]=relColor(d,thr,max)))):TH.cellnull;ctx.fillRect(strip+4+b*cs,a*cs,cs,cs);}}}
  cv.onmousemove=function(e){var rc=cv.getBoundingClientRect(),x=e.clientX-rc.left-strip-4,y=e.clientY-rc.top,a=Math.floor(y/cs),b=Math.floor(x/cs);
    if(a<0||b<0||a>=n||b>=n){tip('');return;}
    var ia=order[a],ib=order[b],sa=r.samples[ia],sb=r.samples[ib],d=relD(r,ia,ib),c=relC(r,ia,ib),g=S.group||{};
    tip('<b>'+esc(sa)+'</b>'+(g[sa]?' ('+esc(g[sa])+')':'')+' &#215; <b>'+esc(sb)+'</b>'+(g[sb]?' ('+esc(g[sb])+')':'')+'<br>'+
      (d==null?'no distance':(d.toLocaleString('en-US')+' SNPs'))+(ia==ib?'':' over '+(c==null?'NA':c+'%')+' of the variable positions'),e.clientX,e.clientY);};
  cv.onmouseleave=function(){tip('');};
  cv.onclick=function(e){var rc=cv.getBoundingClientRect(),a=Math.floor((e.clientY-rc.top)/cs);if(a>=0&&a<n)setHi(r.samples[order[a]]);};
}
function renderRelClusters(){
  var host=el('relclus_body'),S=relSec(),r=relRef(); if(!host)return; if(!(S&&r)){host.innerHTML='';return;}
  var thr=relThr(),minCmp=relMinCmp(),cl=relClusters(r,thr,minCmp).filter(function(c){return c.length>1;}),g=S.group||{},gcol=relGroupColor(),hasG=Object.keys(g).length>0;
  var rows=cl.map(function(c,k){var mx=0,a,b;for(a=0;a<c.length;a++)for(b=a+1;b<c.length;b++){var d=relD(r,c[a],c[b]);if(d!=null&&relOk(r,c[a],c[b],minCmp)&&d>mx)mx=d;}
    var names=c.map(function(i){return r.samples[i];}),gs={};names.forEach(function(s){if(g[s])gs[g[s]]=1;});var ng=Object.keys(gs).length;
    // a large cluster wraps and names its first 40: a 93-sample one made the table 8,000 px wide
    return '<tr><td>'+(k+1)+'</td><td>'+c.length+'</td><td style="text-align:left;white-space:normal" title="'+esc(names.join(', '))+'">'+names.slice(0,40).map(function(s){var col=gcol(s);
        return '<span class="rel-chip"'+(col?' style="border-left:4px solid '+col+'"':'')+' title="'+esc(g[s]||'')+'">'+esc(s)+'</span>';}).join(' ')+(names.length>40?' <span class="na">and '+(names.length-40)+' more</span>':'')+'</td>'+
      '<td>'+mx.toLocaleString('en-US')+'</td>'+(hasG?'<td'+(ng>1?' style="color:var(--warn);font-weight:600" title="the cluster joins samples the samplesheet puts in different groups"':'')+'>'+ng+'</td>':'')+'</tr>';}).join('');
  host.innerHTML='<div class="hot-note" style="margin:10px 16px 6px">Samples joined by a chain of pairs each within '+thr+' SNPs (single linkage). '+
    (hasG?'A cluster spanning several groups is where to look for transmission between them, or for a sample in the wrong group.':'Add a group column to the samplesheet (patient, line, outbreak) to see which clusters cross groups.')+'</div>'+
    '<div class="gtable" style="max-height:48vh"><table><thead><tr><th>#</th><th>Samples</th><th style="text-align:left">Members</th><th>Largest distance</th>'+(hasG?'<th>Groups</th>':'')+'</tr></thead><tbody>'+
    (rows||'<tr><td colspan="5" class="na" style="text-align:left;padding:8px">no two samples are within '+thr+' SNPs.</td></tr>')+'</tbody></table></div>';
  var gh=el('relgroup_body');
  if(gh){var out=R.samples.filter(function(s){return s.grpd;});
    gh.innerHTML=out.length?('<div class="hot-note" style="margin:10px 16px 6px">'+out.length+' '+_plural(out.length,'sample sits','samples sit')+' further than '+(S.threshold)+' SNPs from every other member of '+_plural(out.length,'its group','their groups')+', in groups whose members are otherwise within that of each other. Each is flagged <b>'+esc(flagLab('GROUP_MISMATCH'))+'</b>: a swapped or mislabelled sample, a contaminated or mixed one, or a reinfection. The nearest sample of another group helps tell them apart: a swap sits next to one, a contaminant next to nobody.</div>'+
      '<div class="gtable" style="max-height:40vh"><table><thead><tr><th class="s" style="text-align:left">Sample</th><th style="text-align:left">Its group</th><th>Nearest of its group</th><th>Nearest of another</th></tr></thead><tbody>'+
      out.map(function(s){var d=s.grpd;return '<tr><td class="s" style="text-align:left">'+esc(s.s)+'</td><td style="text-align:left">'+esc(d.g)+'</td>'+
        '<td>'+esc(d.nin)+' &#183; '+d.din.toLocaleString('en-US')+' SNPs</td><td>'+(d.nout?(esc(d.nout)+(d.gout?' ('+esc(d.gout)+')':'')+' &#183; '+d.dout.toLocaleString('en-US')+' SNPs'):'<span class="na">none on its reference</span>')+'</td></tr>';}).join('')+
      '</tbody></table></div>'):'<div class="hot-note" style="margin:10px 16px">'+(Object.keys(S.group||{}).length?'Every sample sits within '+S.threshold+' SNPs of another member of its group, where its group is tight enough to say.':'The samplesheet has no group column (patient, line, series...), so there is nothing to check a sample against.')+'</div>';}
}
