const http = require('http'), fs = require('fs'), path = require('path');
const root = process.argv[2], port = +process.argv[3] || 8765, outFile = process.argv[4];
const inject = fs.readFileSync(path.join(__dirname, 'measure.js'), 'utf8');
const types = {html:'text/html',js:'application/javascript',css:'text/css',png:'image/png',svg:'image/svg+xml',json:'application/json',woff:'font/woff',woff2:'font/woff2'};
http.createServer((req, res) => {
  const url = new URL(req.url, 'http://x');
  if (url.pathname === '/slow.png') { setTimeout(() => { res.writeHead(404); res.end(); }, +(url.searchParams.get('d') || 6000)); return; }
  if (url.pathname === '/report') { let b=''; req.on('data', d => b += d); req.on('end', () => { fs.appendFileSync(outFile, b + '\n'); res.end('ok'); }); return; }
  let p = url.pathname;
  let check = false;
  if (p.endsWith('-check.html')) { p = p.replace('-check.html', '.html'); check = true; }
  const f = path.join(root, decodeURIComponent(p));
  fs.readFile(f, (err, data) => {
    if (err) { res.writeHead(404); res.end(); return; }
    const ext = path.extname(f).slice(1);
    res.writeHead(200, {'Content-Type': types[ext] || 'application/octet-stream'});
    if (check) { const d = url.searchParams.has('measure') ? 40000 : 7000;
      data = data.toString().replace('</body>', '<img src="/slow.png?d=' + d + '" style="display:none"><script>' + inject + '</script></body>'); }
    res.end(data);
  });
}).listen(port, '127.0.0.1');
