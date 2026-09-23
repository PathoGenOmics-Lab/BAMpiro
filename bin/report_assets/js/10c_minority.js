// ---- Minority variants: how low an allele frequency can be trusted ----
// Calls below fixation, the reads they rest on, and, where the samplesheet names libraries of the
// same DNA, how often another library of that DNA calls them too: a real minority is in the DNA,
// an error is not reproduced.
function minBandLabels(M){var e=M.edges;return e.slice(0,-1).map(function(lo,i){return i===0?'&#8804; '+e[1]:lo+'&#8211;'+e[i+1];});}
// The lowest band from which at least 80% of calls reproduce, over bands tested at least 20 times.
function minFloor(M){var r=M.replicates;if(!(r&&r.tested))return null;
  for(var i=0;i<r.tested.length;i++){if(r.tested[i]>=20&&r.reproduced[i]/r.tested[i]>=0.8)return i;}return -1;}
function minBars(vals,labels,fmt,color,title){var mx=Math.max.apply(null,vals.map(function(v){return v||0;}).concat([1e-9]));
  return '<div class="rel-hist">'+vals.map(function(v,i){return '<div class="rel-bar" title="'+esc(title(i))+'"><span class="rel-bv">'+(v==null?'':fmt(v))+'</span>'+
    '<i style="height:'+Math.round(4+(v||0)/mx*70)+'px;background:'+(typeof color=='function'?color(i):color)+'"></i><span class="rel-bl">'+labels[i]+'</span></div>';}).join('')+'</div>';}
function renderMinority(){
  var host=el('minority_body'),M=R.minority; if(!host)return; if(!M){host.innerHTML='';return;}
  var labels=minBandLabels(M),tot=M.hist.reduce(function(a,b){return a+b;},0),n=Object.keys(M.samples).length;
  var html='<div class="hot-note" style="margin:0 0 6px"><b>'+tot.toLocaleString('en-US')+'</b> variant calls below fixation (allele fraction under '+M.fixed+') across '+n+' samples; '+
    '<b>'+(M.with_dp?Math.round(100*M.le3/M.with_dp):0)+'%</b> of them rest on three alternate reads or fewer, where a sequencing error and a real minority look the same.</div>'+
    minBars(M.hist,labels,function(v){return v.toLocaleString('en-US');},'#4f83c2',function(i){return M.hist[i]+' calls at allele fraction '+labels[i].replace('&#8804;','<=').replace('&#8211;','-');});
  var r=M.replicates;
  if(r){
    if(r.tested==null){html+='<div class="hot-note" style="margin:10px 0">The samplesheet names libraries of the same DNA (<code>'+esc(r.column)+'</code>, '+r.pairs+' independent '+_plural(r.pairs,'pair')+'), '+
      'but without the SNP matrix a call another library missed cannot be told from a site it did not read, so they are not compared.</div>';}
    else{
      var rates=r.tested.map(function(t,i){return t?r.reproduced[i]/t:null;}),fl=minFloor(M),fx=r.fixed_tested?r.fixed_reproduced/r.fixed_tested:null;
      html+='<div class="dsub" style="margin:14px 0 4px">Reproduced by another library of the same DNA <span style="font-weight:400;color:#94a3b8">('+
        esc(r.column)+'; '+r.pairs+' independent '+_plural(r.pairs,'pair')+(r.shared?(', '+r.shared+' left out for sharing reads'):'')+')</span></div>'+
        minBars(rates,labels,function(v){return Math.round(v*100)+'%';},function(i){return rates[i]!=null&&rates[i]>=0.8?'#22a06b':(rates[i]!=null&&rates[i]>=0.5?'#e6b25a':'#e0544f');},
          function(i){return r.reproduced[i]+' of '+r.tested[i]+' calls reproduced';})+
        '<div class="hot-note">'+(fl==null||fl<0?'No band of allele fraction reaches 80% reproduced'+(fx!=null?'; fixed calls reach '+Math.round(fx*100)+'%':'')+'.':
          ((fl===0?'At every allele fraction, 80% or more of the calls are reproduced':'From allele fraction '+M.edges[fl]+' up, 80% or more of the calls are reproduced')+(fl>0?('; below it, '+Math.round(100*(rates.slice(0,fl).reduce(function(a,v,i){return a+(v||0)*r.tested[i];},0)/Math.max(1,r.tested.slice(0,fl).reduce(function(a,b){return a+b;},0))))+'% are'):'')+
          (fx!=null?'. Fixed calls: '+Math.round(fx*100)+'%, the ceiling':'')+'.'))+
        ' Only calls the other library was read at count; a site it did not read says nothing.</div>';}}
  else html+='<div class="hot-note" style="margin:10px 0">Add a column naming each sample&#39;s DNA extract to the samplesheet (<code>dna_id</code>, <code>extract</code>, <code>biosample</code>...) and libraries of the same DNA are compared to measure where the noise ends.</div>';
  var rows=Object.keys(M.samples).map(function(s){var x=M.samples[s];return {s:s,x:x};}).sort(function(a,b){return b.x.n-a.x.n;});
  var big=el('minorityPanel')&&el('minorityPanel').classList.contains('expanded'),cap=big?rows.length:25;
  html+='<div class="gtable" style="max-height:'+(big?'60vh':'38vh')+';margin-top:10px"><table><thead><tr><th class="s" style="text-align:left">Sample</th><th>Calls below fixation</th><th>On 3 reads or fewer</th><th>Median alternate reads</th>'+
    labels.map(function(l){return '<th>'+l+'</th>';}).join('')+'</tr></thead><tbody>'+
    rows.slice(0,cap).map(function(o){var x=o.x;return '<tr><td class="s" style="text-align:left">'+esc(o.s)+'</td><td>'+x.n.toLocaleString('en-US')+'</td>'+
      '<td>'+(x.with_dp?Math.round(100*x.le3/x.with_dp)+'%':'NA')+'</td><td>'+(x.median_alt==null?'NA':x.median_alt)+'</td>'+
      x.hist.map(function(c){return '<td>'+(c||'')+'</td>';}).join('')+'</tr>';}).join('')+'</tbody></table></div>'+
    (rows.length>cap?'<div class="hot-note">Showing '+cap+' of '+rows.length+' samples; open the panel full for the rest.</div>':'');
  host.innerHTML=html;
}
