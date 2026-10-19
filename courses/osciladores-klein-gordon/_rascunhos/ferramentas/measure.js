(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  for (let i = 0; i < 100; i++) {
    if (window.MathJax && MathJax.Hub) { await new Promise(r => MathJax.Hub.Queue(r)); break; }
    await wait(100);
  }
  await wait(1000);
  const ms = location.search.match(/slide=(\d+)/);
  if (ms && window.Reveal) {
    const before = { ready: Reveal.isReady(), idx: Reveal.getIndices(), hash: location.hash, t: performance.now() };
    if (!Reveal.isReady()) { await new Promise(r => Reveal.on('ready', r)); }
    const ix = Reveal.getIndices(Reveal.getSlides()[+ms[1]]);
    Reveal.slide(ix.h, ix.v);
    await wait(800);
    const after = { ready: Reveal.isReady(), idx: Reveal.getIndices(), hash: location.hash, t: performance.now() };
    await fetch('/report', { method: 'POST', body: JSON.stringify({ before, after }) });
    return;
  }
  if (!window.Reveal || !location.search.includes('measure')) return;
  // reveal.js lazy-loads data-src images near the current slide only; load them all first.
  await Promise.all(Array.from(document.querySelectorAll('.reveal .slides img[data-src]')).map(img => new Promise(r => {
    if (img.getAttribute('src')) return r();
    img.onload = img.onerror = () => r();
    img.setAttribute('src', img.getAttribute('data-src'));
  })));
  const out = [];
  const scale = Reveal.getScale();
  const secs = Array.from(document.querySelectorAll('.reveal .slides section')).filter(s => !s.querySelector('section'));
  secs.forEach((sec, i) => {
    const chain = [sec];
    if (sec.parentElement && sec.parentElement.tagName === 'SECTION') chain.push(sec.parentElement);
    const prevs = chain.map(el => ({ el, d: el.style.display, v: el.style.visibility, o: el.style.opacity }));
    chain.forEach(el => { el.style.display = 'block'; el.style.visibility = 'visible'; el.style.opacity = '1'; });
    const r = sec.getBoundingClientRect();
    let maxRight = r.right, maxBottom = r.bottom, worst = '', worstB = '';
    sec.querySelectorAll('*').forEach(el => {
      const b = el.getBoundingClientRect();
      if (b.width === 0 || el.closest('aside') || el.closest('.MJX_Assistive_MathML')) return;
      if (b.right > maxRight) { maxRight = b.right; worst = el.tagName + (el.className ? '.' + String(el.className).split(' ')[0] : ''); }
      if (b.bottom > maxBottom) { maxBottom = b.bottom; worstB = el.tagName + (el.className ? '.' + String(el.className).split(' ')[0] : '') + (el.id ? '#' + el.id : ''); }
    });
    const h = sec.querySelector('h1,h2,h3');
    const cfg = Reveal.getConfig();
    const usable = cfg.height * (1 - 2 * cfg.margin);
    const contentH = (maxBottom - r.top) / scale;
    out.push({ i, title: h ? h.textContent.trim().slice(0, 40) : '', secW: Math.round(r.width / scale), secH: Math.round(contentH), overflowRight: Math.round((maxRight - r.right) / scale), overflowBottom: Math.round(Math.max(0, contentH - usable)), worst, worstB, scale: +scale.toFixed(3) });
    prevs.forEach(p => { p.el.style.display = p.d; p.el.style.visibility = p.v; p.el.style.opacity = p.o; });
  });
  await fetch('/report', { method: 'POST', body: JSON.stringify(out) });
})();
