#!/usr/bin/env node
/**
 * snapshot-driver.mjs — THE ONE DRIVER of the Cortical Flow viewer's snapshot panel, for rheome.snapshot and nxr.snapshot (Rheome's
 * `+rheome/private/snapshot_driver.mjs`, generalised; shipped in the viewer's build at `driver/snapshot-driver.mjs`). Headless Chrome
 * over raw DevTools (no dependencies: Node 22's fetch and WebSocket), the viewer and one bundle served on 127.0.0.1.
 *
 *   node snapshot-driver.mjs job.json
 *
 * A SNAPSHOT JOB (`nxr-snapshot-job@1`; relative paths are relative to the job file):
 *   { "format": "nxr-snapshot-job@1",
 *     "viewer": <dist-viewer folder>   (default: the build this driver sits in, `..`),
 *     "bundle": <rheome.bundle@1 folder>                                                     REQUIRED
 *     "out":    <folder written>                                                              REQUIRED
 *     "chrome": <Chrome or Chromium>   (default $RHEOME_CHROME, $CHROME, else the usual install),
 *     "tokens": <figure-tokens.json served as figure/figure-tokens.json>, "tokens_sha256", "tokens_path" (written into each marks file),
 *     "size": [w, h] (the window, default [1000, 760]), "wait": ms after the bundle opens (default 600), "timeout": s per shot (default 120),
 *     "shots": [{ "name": "panel", "shot": <the viewer's Shot: src/viewer/print/panel.ts> }, …] }                       REQUIRED
 *   → for each shot: <out>/<name>.svg (the panel), <name>.marks.json (nxr.marks@1), <name>.report.json (rheome.snapshot's `out`, 0-based;
 *     `build` the release's VERSION when the viewer folder has one), and each view's raw render <out>/<name>/<label>.png; <out>/driver.json
 *     the viewer's `ready` and each shot's warnings.
 * A PRINT JOB, as Rheome's driver took it: { "chrome", "print": [{ svg, pdf, png, width_mm, height_mm, dpi }] }: each panel printed to PDF
 *   (vector text, exact page size) and PNG (at dpi).
 *
 * Every reply the page gives (the SVG, each PNG, the report, the marks) stays in the page and is read back in 1 MB slices: a DevTools
 * reply over ~4 MB never arrives. Exit status 1, with the viewer's own words, on anything refused; nothing is half written for a shot.
 */
import { spawn } from 'node:child_process'
import { createHash } from 'node:crypto'
import fs from 'node:fs'
import http from 'node:http'
import os from 'node:os'
import path from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'

export const FORMAT = 'nxr-snapshot-job@1'
export const SLICE = 1 << 20
const HERE = path.dirname(fileURLToPath(import.meta.url))

/** The job, checked before a browser starts; relative paths made absolute (against the job file's folder). */
export function parseJob(text, jobPath = 'job.json') {
  let j
  try { j = JSON.parse(text) } catch (e) { throw new Error(`${jobPath}: not valid JSON (${e.message})`) }
  if (!j || typeof j !== 'object' || Array.isArray(j)) throw new Error(`${jobPath}: a job is a JSON object`)
  const dir = path.dirname(path.resolve(jobPath)), abs = (p) => (typeof p === 'string' && p ? path.resolve(dir, p) : p)
  if (j.print !== undefined) {
    if (!Array.isArray(j.print) || !j.print.every((p) => p && typeof p.svg === 'string' && p.width_mm > 0 && p.height_mm > 0))
      throw new Error(`${jobPath}: "print" is a list of {svg, pdf?, png?, width_mm, height_mm, dpi?}`)
    return { ...j, chrome: chromePath(j.chrome), print: j.print.map((p) => ({ ...p, svg: abs(p.svg), pdf: abs(p.pdf), png: abs(p.png) })) }
  }
  if (j.format !== FORMAT) throw new Error(`${jobPath}: "format" must be "${FORMAT}"`)
  for (const k of ['bundle', 'out']) if (typeof j[k] !== 'string' || !j[k]) throw new Error(`${jobPath}: "${k}" is required (a path)`)
  if (!Array.isArray(j.shots) || !j.shots.length) throw new Error(`${jobPath}: "shots" is a list of {name, shot}`)
  const names = new Set()
  for (const [i, s] of j.shots.entries()) {
    if (!s || typeof s !== 'object' || typeof s.name !== 'string' || !/^[\w\-.]+$/.test(s.name)) throw new Error(`${jobPath}: shots[${i}].name: letters, digits, - _ . only`)
    if (names.has(s.name)) throw new Error(`${jobPath}: shots[${i}].name "${s.name}" is used twice`)
    names.add(s.name)
    if (!s.shot || typeof s.shot !== 'object' || Array.isArray(s.shot)) throw new Error(`${jobPath}: shots[${i}].shot is the viewer's Shot (an object)`)
    if (typeof s.shot.field !== 'string' || typeof s.shot.unit !== 'string' || !(s.shot.widthMm > 0)) throw new Error(`${jobPath}: shots[${i}].shot needs field, unit and widthMm`)
  }
  if (j.size !== undefined && !(Array.isArray(j.size) && j.size.length === 2 && j.size.every((n) => Number.isInteger(n) && n > 0))) throw new Error(`${jobPath}: "size" is [width, height] in pixels`)
  return { ...j, viewer: abs(j.viewer) ?? path.resolve(HERE, '..'), bundle: abs(j.bundle), out: abs(j.out), tokens: abs(j.tokens), chrome: chromePath(j.chrome) }
}

/** Chrome: given, $RHEOME_CHROME, $CHROME, else the usual install */
export function chromePath(c) {
  if (c) return c
  if (process.env.RHEOME_CHROME) return process.env.RHEOME_CHROME
  if (process.env.CHROME) return process.env.CHROME
  const usual = { darwin: ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '/Applications/Chromium.app/Contents/MacOS/Chromium'],
    linux: ['/usr/bin/google-chrome', '/usr/bin/chromium', '/usr/bin/chromium-browser'], win32: ['C:/Program Files/Google/Chrome/Application/chrome.exe'] }[process.platform] ?? []
  return usual.find((p) => fs.existsSync(p)) ?? usual[0] ?? 'chrome'
}

/** A string held by the page (`expr` evaluates to it there), read back in `slice`-character pieces through `evalp` (one DevTools
 *  evaluation each): a reply over ~4 MB never arrives whole. */
export async function readSliced(evalp, expr, slice = SLICE) {
  const n = await evalp(`(${expr}).length`)
  if (!Number.isInteger(n) || n < 0) throw new Error(`the page holds no string at ${expr}`)
  const parts = []
  for (let i = 0; i < n; i += slice) parts.push(await evalp(`(${expr}).slice(${i}, ${i + slice})`))
  const s = parts.join('')
  if (s.length !== n) throw new Error(`read ${s.length} of ${n} characters of ${expr}`)
  return s
}

/** the release's build (a VERSION file beside viewer.html: `build: …`), or null for a development build */
export function releaseBuild(viewer) {
  const f = path.join(viewer, 'VERSION')
  if (!fs.existsSync(f)) return null
  return /^build:\s*(\S+)/m.exec(fs.readFileSync(f, 'utf8'))?.[1] ?? null
}

const sleep = (ms) => new Promise((r) => setTimeout(r, ms))
let chrome, server, prof
function cleanup() {
  try { chrome?.kill('SIGKILL') } catch {}
  try { server?.close() } catch {}
  try { if (prof) fs.rmSync(prof, { recursive: true, force: true }) } catch {}
}
function fail(msg) { console.error('snapshot: ' + msg); cleanup(); process.exit(1) }

async function devtools(job, W, H) {
  prof = fs.mkdtempSync(path.join(os.tmpdir(), 'nxr-snapshot-chrome-'))
  chrome = spawn(job.chrome, ['--headless=new', '--remote-debugging-port=0', `--user-data-dir=${prof}`, '--enable-unsafe-webgpu',
    ...(process.platform === 'darwin' ? ['--use-angle=metal'] : []), '--ignore-gpu-blocklist',
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

/** each composed panel printed to PDF and PNG (Rheome's print job, unchanged) */
async function printPanels(job) {
  const { ws, send, evalp } = await devtools(job, 800, 600)
  setTimeout(() => fail('printing timed out'), 1000 * (job.timeout ?? 120))
  for (const p of job.print) {
    const wpx = Math.ceil(p.width_mm / 25.4 * 96), hpx = Math.ceil(p.height_mm / 25.4 * 96)
    await send('Emulation.clearDeviceMetricsOverride')
    const html = p.svg.replace(/\.svg$/, '') + '.print.html'             // inline in an HTML page (printToPDF does not finish on a bare SVG document)
    fs.writeFileSync(html, `<!doctype html><html><head><meta charset="utf-8"><style>@page{size:${p.width_mm}mm ${p.height_mm}mm;margin:0}html,body{margin:0;padding:0}svg{display:block}</style></head><body>${fs.readFileSync(p.svg, 'utf8').replace(/^<\?xml[^>]*>\s*/, '')}</body></html>`)
    await send('Page.navigate', { url: pathToFileURL(html).href })
    for (let ok = false; !ok;) { await sleep(100); ok = await evalp(`document.readyState === 'complete' && location.href.startsWith('file:')`) }
    await sleep(200)
    if (p.pdf) {
      const r = await send('Page.printToPDF', { paperWidth: p.width_mm / 25.4, paperHeight: p.height_mm / 25.4, marginTop: 0, marginBottom: 0,
        marginLeft: 0, marginRight: 0, printBackground: true, pageRanges: '1', scale: 1, preferCSSPageSize: true, transferMode: 'ReturnAsStream' })
      if (!r.result?.stream) fail(`${p.svg}: no PDF (${JSON.stringify(r.error ?? {})})`)
      const parts = []                                                     // in 1 MB reads: a PDF of several MB in one reply does not arrive
      for (let eof = false; !eof;) {
        const c = (await send('IO.read', { handle: r.result.stream, size: SLICE })).result
        parts.push(Buffer.from(c.data, c.base64Encoded ? 'base64' : 'utf8')); eof = c.eof
      }
      await send('IO.close', { handle: r.result.stream })
      fs.writeFileSync(p.pdf, Buffer.concat(parts))
    }
    if (p.png) {                                                           // the PNG at dpi: the page drawn at dpi / 96 device pixels per CSS pixel
      await send('Emulation.setDeviceMetricsOverride', { width: wpx, height: hpx, deviceScaleFactor: (p.dpi ?? 600) / 96, mobile: false })
      await sleep(200)
      const r = await send('Page.captureScreenshot', { format: 'png', clip: { x: 0, y: 0, width: p.width_mm / 25.4 * 96, height: p.height_mm / 25.4 * 96, scale: 1 } })
      if (!r.result?.data) fail(`${p.svg}: no PNG`)
      fs.writeFileSync(p.png, Buffer.from(r.result.data, 'base64'))
    }
    fs.rmSync(html, { force: true })
  }
  ws.close(); cleanup(); process.exit(0)
}

const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.json': 'application/json',
  '.wasm': 'application/wasm', '.svg': 'image/svg+xml', '.png': 'image/png', '.woff2': 'font/woff2' }

/** the snapshot shots: the viewer opens the bundle; each shot is the viewer's `figure(shot)`, its answer read back in slices */
async function snapshots(job) {
  if (!fs.existsSync(path.join(job.viewer, 'viewer.html'))) fail(`no viewer at ${job.viewer} (viewer.html)`)
  if (!fs.existsSync(path.join(job.bundle, 'manifest.json'))) fail(`no bundle at ${job.bundle} (manifest.json)`)
  const [W, H] = job.size ?? [1000, 760]
  /* the files: /b/<file> from the bundle, everything else from the viewer (read only); the tokens the job names, not the viewer's copy */
  server = http.createServer((req, res) => {
    const u = decodeURIComponent(new URL(req.url, 'http://x').pathname)
    if (job.tokens && u === '/figure/figure-tokens.json') {
      if (!fs.existsSync(job.tokens)) { res.writeHead(500); res.end(); fail('the figure tokens are missing: ' + job.tokens); return }
      res.writeHead(200, { 'content-type': 'application/json' }); res.end(fs.readFileSync(job.tokens)); return
    }
    const [root, rel] = u.startsWith('/b/') ? [job.bundle, u.slice(3)] : [job.viewer, u.slice(1) || 'viewer.html']
    const f = path.resolve(root, rel)
    if (!f.startsWith(path.resolve(root) + path.sep) || !fs.existsSync(f) || !fs.statSync(f).isFile()) { res.writeHead(404); res.end(); return }
    res.writeHead(200, { 'content-type': TYPES[path.extname(f)] ?? 'application/octet-stream' })
    fs.createReadStream(f).pipe(res)
  })
  await new Promise((r) => server.listen(0, '127.0.0.1', r))
  const url = `http://127.0.0.1:${server.address().port}/viewer.html?bundle=${encodeURIComponent('/b/manifest.json')}`
  const { ws, send, evalp, logs } = await devtools(job, W, H)
  await send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: 1, mobile: false })
  await send('Page.addScriptToEvaluateOnNewDocument', { source:
    `window.__ev=[];for(const n of ['ready','error','warning','figure'])window.addEventListener('viewer:'+n,e=>window.__ev.push([n,e.detail]))` })
  await send('Page.navigate', { url })
  let opened = null
  for (let t = Date.now(); !opened && Date.now() - t < 1000 * (job.timeout ?? 120);) {
    await sleep(100)
    opened = await evalp(`(window.__ev||[]).find(e=>e[0]==='ready'||(e[0]==='error'&&e[1].code!=='settings'))||null`)
  }
  if (!opened) fail('the viewer did not open the bundle in time\n' + logs.join('\n'))
  if (opened[0] === 'error') fail(`the viewer: ${opened[1].code}: ${opened[1].message}`)
  await sleep(job.wait ?? 600)
  fs.mkdirSync(job.out, { recursive: true })
  const build = releaseBuild(job.viewer), done = { ready: opened[1], build, shots: [] }
  for (const { name, shot } of job.shots) {
    let timer
    /* the figure stays in the page (window.__fig: svg, marks and report as JSON text, each render's PNG); only its sizes come back here */
    const head = await Promise.race([evalp(`(async()=>{window.__ev=[];const p=new Promise(r=>{window.addEventListener('viewer:figure',e=>r({ok:e.detail}),{once:true});`
      + `window.addEventListener('viewer:error',e=>r({err:e.detail}),{once:true})});await window.viewer.figure(${JSON.stringify(shot)}, true);const a=await p;`
      + `if(a.err)return{error:a.err.message};const d=a.ok;window.__fig={svg:d.svg,marks:JSON.stringify(d.marks),report:JSON.stringify(d.report),renders:d.renders||[]};`
      + `return{renders:window.__fig.renders.map(r=>r.label),warnings:window.__ev.filter(e=>e[0]==='warning').map(e=>e[1].message)}})()`),
    new Promise((r) => { timer = setTimeout(r, 1000 * (job.timeout ?? 120)) }).then(() => fail(`${name}: no figure in ${job.timeout ?? 120} s\n${logs.slice(-20).join('\n')}`))])
    clearTimeout(timer)
    if (!head || head.error) fail(`${name}: the viewer refused the shot: ${head?.error ?? 'no answer'}`)
    const svg = await readSliced(evalp, 'window.__fig.svg'), marks = JSON.parse(await readSliced(evalp, 'window.__fig.marks'))
    const report = JSON.parse(await readSliced(evalp, 'window.__fig.report'))
    const pngs = []
    for (const [i, label] of head.renders.entries()) {
      const png = await readSliced(evalp, `window.__fig.renders[${i}].png`)
      if (!png.startsWith('data:image/png;base64,')) fail(`${name}: the render ${label} is no PNG`)
      pngs.push([label, Buffer.from(png.slice(png.indexOf(',') + 1), 'base64')])
    }
    await evalp('delete window.__fig')
    if (build) report.build = build                                        // a release records its VERSION's build
    const dir = path.join(job.out, name)
    fs.mkdirSync(dir, { recursive: true })
    for (const [label, buf] of pngs) fs.writeFileSync(path.join(dir, `${label}.png`), buf)
    fs.writeFileSync(path.join(job.out, `${name}.svg`), svg)
    const sha = job.tokens_sha256 ?? (job.tokens && fs.existsSync(job.tokens) ? createHash('sha256').update(fs.readFileSync(job.tokens)).digest('hex') : undefined)
    fs.writeFileSync(path.join(job.out, `${name}.marks.json`), JSON.stringify(job.tokens ? { ...marks, tokens_sha256: sha, tokens_path: job.tokens_path ?? job.tokens } : marks, null, 2))
    fs.writeFileSync(path.join(job.out, `${name}.report.json`), JSON.stringify(report, null, 2))
    done.shots.push({ name, renders: head.renders, warnings: head.warnings })
    for (const w of head.warnings) console.error(`snapshot: ${name}: ${w}`)
  }
  fs.writeFileSync(path.join(job.out, 'driver.json'), JSON.stringify(done, null, 2))
  ws.close(); cleanup(); process.exit(0)
}

if (process.argv[1] && import.meta.url === pathToFileURL(path.resolve(process.argv[1])).href) {
  process.on('exit', cleanup)
  const file = process.argv[2]
  if (!file) fail('usage: node snapshot-driver.mjs <job.json>')
  let job
  try { job = parseJob(fs.readFileSync(file, 'utf8'), file) } catch (e) { fail(e.message) }
  if (job.print) await printPanels(job)
  else await snapshots(job)
}
