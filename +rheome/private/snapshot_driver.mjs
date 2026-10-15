// The Cortical Flow viewer (dist-viewer) in headless Chrome, for rheome.snapshot: serve the viewer and one bundle, open it, apply
// each shot's settings, save the viewer's own screenshot (PNG, alpha kept) per shot. Graphics' F1-board shoot.mjs, made a tool.
// Usage: node snapshot_driver.mjs <job.json>
// job: { viewer: <dist-viewer folder>, bundle: <bundle folder>, out: <folder>, chrome: <path>, size: [w,h], scale: 2, wait: 600,
//        timeout: 60, tokens: <figure-tokens.json to serve as figure/figure-tokens.json>, tokens_sha256, shots: [{ name, settings, scale? }] }    -> <out>/<name>.png each, and the events in <out>/driver.json
//   or { chrome, print: [{ svg, pdf, png, width_mm, height_mm, dpi }] }: each composed panel (an SVG, views embedded at their
//        resolution) printed to PDF (vector text, exact page size) and PNG (at dpi). MATLAB's own PDF export resamples images
//        to screen pixels, which is why the panel is printed here.
import { spawn } from 'node:child_process'
import fs from 'node:fs'
import http from 'node:http'
import os from 'node:os'
import path from 'node:path'

const job = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'))
const [W, H] = job.size ?? [1000, 760]
const sleep = (ms) => new Promise((r) => setTimeout(r, ms))
const fail = (msg) => { console.error('snapshot: ' + msg); cleanup(); process.exit(1) }
let chrome, server, prof
function cleanup() {
  try { chrome?.kill('SIGKILL') } catch {}
  try { server?.close() } catch {}
  try { if (prof) fs.rmSync(prof, { recursive: true, force: true }) } catch {}
}
process.on('exit', cleanup)

async function devtools(W, H) {
  prof = fs.mkdtempSync(path.join(os.tmpdir(), 'rheome-chrome-'))
  chrome = spawn(job.chrome, ['--headless=new', '--remote-debugging-port=0', `--user-data-dir=${prof}`, '--enable-unsafe-webgpu',
    '--no-first-run', '--no-default-browser-check', `--window-size=${W},${H}`, '--hide-scrollbars', '--allow-file-access-from-files', 'about:blank'],
    { stdio: ['ignore', 'ignore', 'pipe'] })
  chrome.on('error', (e) => fail(`cannot start Chrome (${job.chrome}): ${e.message}`))
  const wsUrl = await new Promise((ok) => {
    let buf = ''
    chrome.stderr.on('data', (d) => { buf += d; const m = buf.match(/DevTools listening on (ws:\/\/\S+)/); if (m) ok(m[1]) })
    setTimeout(() => ok(null), 20000)
  })
  if (!wsUrl) fail('Chrome did not open its DevTools port')
  const port = new URL(wsUrl).port
  const page = (await (await fetch(`http://127.0.0.1:${port}/json/list`)).json()).find((t) => t.type === 'page')
  const ws = new WebSocket(page.webSocketDebuggerUrl)
  await new Promise((r, j) => { ws.addEventListener('open', r); ws.addEventListener('error', j) })
  let id = 0; const pend = new Map(); const logs = []
  ws.addEventListener('message', (ev) => {
    const m = JSON.parse(ev.data)
    if (m.id && pend.has(m.id)) { pend.get(m.id)(m); pend.delete(m.id) }
    else if (m.method === 'Runtime.consoleAPICalled') logs.push(m.params.args.map((a) => a.value ?? a.description).join(' '))
    else if (m.method === 'Runtime.exceptionThrown') logs.push('EXC ' + JSON.stringify(m.params.exceptionDetails).slice(0, 400))
  })
  const send = (method, params = {}) => new Promise((r) => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })) })
  const evalp = async (expr) => {
    const r = await send('Runtime.evaluate', { expression: expr, awaitPromise: true, returnByValue: true })
    if (r.result?.exceptionDetails) throw new Error(JSON.stringify(r.result.exceptionDetails).slice(0, 600))
    return r.result?.result?.value
  }
  await send('Runtime.enable'); await send('Page.enable')
  return { ws, send, evalp, logs }
}

async function printPanels() {
  const { ws, send, evalp } = await devtools(800, 600)
  setTimeout(() => fail('printing timed out'), 1000 * (job.timeout ?? 120))
  for (const p of job.print) {
    const wpx = Math.ceil(p.width_mm / 25.4 * 96), hpx = Math.ceil(p.height_mm / 25.4 * 96)
    await send('Emulation.clearDeviceMetricsOverride')
    // inline in an HTML page (printToPDF does not finish on a bare SVG document)
    const html = p.svg.replace(/\.svg$/, '') + '.print.html'
    fs.writeFileSync(html, `<!doctype html><html><head><meta charset="utf-8"><style>@page{size:${p.width_mm}mm ${p.height_mm}mm;margin:0}html,body{margin:0;padding:0}svg{display:block}</style></head><body>${fs.readFileSync(p.svg, 'utf8')}</body></html>`)
    await send('Page.navigate', { url: 'file://' + html })
    for (let ok = false; !ok;) { await sleep(100); ok = await evalp(`document.readyState === 'complete' && location.href.startsWith('file:')`) }
    await sleep(200)
    if (p.pdf) {
      const r = await send('Page.printToPDF', { paperWidth: p.width_mm / 25.4, paperHeight: p.height_mm / 25.4, marginTop: 0, marginBottom: 0,
        marginLeft: 0, marginRight: 0, printBackground: true, pageRanges: '1', scale: 1, preferCSSPageSize: true, transferMode: 'ReturnAsStream' })
      if (!r.result?.stream) fail(`${p.svg}: no PDF (${JSON.stringify(r.error ?? {})})`)
      const parts = []   // read in chunks: a PDF of several MB in one reply does not arrive
      for (let eof = false; !eof;) {
        const c = (await send('IO.read', { handle: r.result.stream, size: 1 << 20 })).result
        parts.push(Buffer.from(c.data, c.base64Encoded ? 'base64' : 'utf8')); eof = c.eof
      }
      await send('IO.close', { handle: r.result.stream })
      fs.writeFileSync(p.pdf, Buffer.concat(parts))
    }
    if (p.png) {   // the PNG at dpi: the page drawn at dpi / 96 device pixels per CSS pixel
      await send('Emulation.setDeviceMetricsOverride', { width: wpx, height: hpx, deviceScaleFactor: p.dpi / 96, mobile: false })
      await sleep(200)
      const r = await send('Page.captureScreenshot', { format: 'png', clip: { x: 0, y: 0, width: p.width_mm / 25.4 * 96, height: p.height_mm / 25.4 * 96, scale: 1 } })
      if (!r.result?.data) fail(`${p.svg}: no PNG`)
      fs.writeFileSync(p.png, Buffer.from(r.result.data, 'base64'))
    }
    fs.rmSync(html, { force: true })
  }
  ws.close(); cleanup(); process.exit(0)
}

if (job.print) await printPanels()
/* the files: /b/<file> from the bundle, everything else from the viewer (read only); .bin as octets (the viewer's types) */
const types = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.json': 'application/json',
  '.wasm': 'application/wasm', '.svg': 'image/svg+xml', '.png': 'image/png', '.woff2': 'font/woff2' }
server = http.createServer((req, res) => {
  const u = decodeURIComponent(new URL(req.url, 'http://x').pathname)
  if (job.tokens && u === '/figure/figure-tokens.json') {   // the viewer's fetch("./figure/figure-tokens.json")
    if (!fs.existsSync(job.tokens)) { res.writeHead(500); res.end(); fail('the figure tokens are missing: ' + job.tokens); return }
 res.writeHead(200, { 'content-type': 'application/json' }); res.end(fs.readFileSync(job.tokens)); return }   // the tokens used, not the viewer's copy
  const [root, rel] = u.startsWith('/b/') ? [job.bundle, u.slice(3)] : [job.viewer, u.slice(1) || 'viewer.html']
  const f = path.resolve(root, rel)
  if (!f.startsWith(path.resolve(root) + path.sep) || !fs.existsSync(f) || !fs.statSync(f).isFile()) { res.writeHead(404); res.end(); return }
  res.writeHead(200, { 'content-type': types[path.extname(f)] ?? 'application/octet-stream' })
  fs.createReadStream(f).pipe(res)
})
await new Promise((r) => server.listen(0, '127.0.0.1', r))
const url = `http://127.0.0.1:${server.address().port}/viewer.html?bundle=${encodeURIComponent('/b/manifest.json')}`

const { ws, send, evalp, logs } = await devtools(W, H)
await send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: 1, mobile: false })
await send('Page.addScriptToEvaluateOnNewDocument', { source:
  `window.__ev=[];for(const n of ['ready','error','screenshot','settings'])window.addEventListener('viewer:'+n,e=>window.__ev.push([n,e.detail]))` })
await send('Page.navigate', { url })
let opened = null
for (let t = Date.now(); !opened && Date.now() - t < 1000 * (job.timeout ?? 60);) {
  await sleep(100)
  opened = await evalp(`(window.__ev||[]).find(e=>e[0]==='ready'||(e[0]==='error'&&e[1].code!=='settings'))`)
}
if (!opened) fail('the viewer did not open the bundle in time\n' + logs.join('\n'))
if (opened[0] === 'error') fail(`the viewer: ${opened[1].code}: ${opened[1].message}`)

const done = { ready: opened[1], shots: [] }
let timer
for (const shot of job.shots) {
  await evalp(`(async()=>{window.__ev=[];window.viewer.set(${JSON.stringify(shot.settings)})})()`)
  await sleep(shot.wait ?? job.wait ?? 600)
  const errs = await evalp(`window.__ev.filter(e=>e[0]==='error').map(e=>e[1])`)
  if (errs?.length) fail(`${shot.name}: the viewer refused ${JSON.stringify(errs)}`)
  // a screenshot the viewer never answers (e.g. a canvas too large for the GPU) fails the run instead of hanging it
  const shotEv = await Promise.race([evalp(`(async()=>{const p=new Promise(r=>{window.addEventListener('viewer:screenshot',e=>r({png:e.detail.png,marks:e.detail.marks}),{once:true});` +
    `window.addEventListener('viewer:error',e=>r({png:'error: '+e.detail.code+': '+e.detail.message}),{once:true})});` +
    `await window.viewer.screenshot(${shot.scale ?? job.scale ?? 2}, true);return await p})()`),
    new Promise((r) => { timer = setTimeout(r, 1000 * (job.timeout ?? 60)) })
      .then(() => fail(`${shot.name}: no screenshot in ${job.timeout ?? 60} s at pixel ratio ${shot.scale ?? job.scale ?? 2}`))])
  clearTimeout(timer)   // this shot's clock only
  const png = shotEv?.png
  if (!png?.startsWith('data:image/png;base64,')) fail(`${shot.name}: no screenshot${png ? ' (' + png + ')' : ''}`)
  fs.writeFileSync(path.join(job.out, shot.name + '.png'), Buffer.from(png.split(',')[1], 'base64'))
  if (shotEv.marks) fs.writeFileSync(path.join(job.out, shot.name + '.marks.json'), JSON.stringify(job.tokens ? { ...shotEv.marks, tokens_sha256: job.tokens_sha256, tokens_path: job.tokens_path ?? job.tokens } : shotEv.marks, null, 2))   // keyed overlays
  done.shots.push({ name: shot.name, settings: await evalp('window.viewer.settings') })
}
fs.writeFileSync(path.join(job.out, 'driver.json'), JSON.stringify(done))
ws.close(); cleanup(); process.exit(0)
