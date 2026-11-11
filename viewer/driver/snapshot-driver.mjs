#!/usr/bin/env node
/**
 * snapshot-driver.mjs — THE ONE DRIVER of the snapshot panel, for rheome.snapshot and nxr.snapshot (Rheome's
 * `+rheome/private/snapshot_driver.mjs`, generalised; shipped in the Editor's build at `driver/snapshot-driver.mjs`). Headless Chrome
 * over raw DevTools (no dependencies: Node 22's fetch and WebSocket), the Editor and one bundle served on 127.0.0.1; the Editor draws
 * the panels (editor/figure/cortexPanel.ts). The viewer (`viewer.html`) that drew them before was retired on 2026-11-11.
 *
 *   node snapshot-driver.mjs job.json
 *
 * A SNAPSHOT JOB (`nxr-snapshot-job@1`; relative paths are relative to the job file):
 *   { "format": "nxr-snapshot-job@1",
 *     "editor": <the Editor's build folder, or a running Editor's url>   (default: the build this driver sits in, `..`),
 *     "bundle": <rheome.bundle@1 folder>                                                     REQUIRED
 *     "out":    <folder written>                                                              REQUIRED
 *     "chrome": <Chrome or Chromium>   (default $RHEOME_CHROME, $CHROME, else the usual install),
 *     "tokens": <figure-tokens.json served as figure/figure-tokens.json>, "tokens_sha256", "tokens_path" (written into each marks file),
 *     "size": [w, h] (the window, default [1000, 760]), "wait": ms after the bundle opens (default 600), "timeout": s per shot (default 120),
 *     "shots": [{ "name": "panel", "shot": <the Shot: editor/src/editor/figure/print/panel.ts> }, …] }                       REQUIRED
 *   → for each shot: <out>/<name>.svg (the panel), <name>.marks.json (nxr.marks@1), <name>.report.json (rheome.snapshot's `out`, 0-based;
 *     `build` the release's VERSION when the Editor's build has one), and each view's raw render <out>/<name>/<label>.png; <out>/driver.json
 *     the viewer's `ready` and each shot's warnings.
 *   RAW PICTURES (a shot with `stagePx` [w, h], panel.ts `Shot`): each render exactly w × h, and beside it <out>/<name>/<label>.marks.json
 *     (nxr.marks@1: that render's overlays, pixels in its raster); the report adds `stops`, `fit` and `renders` (each render's view, frame,
 *     size, mm per pixel, fit, colourbar); <name>.svg only when the shot gives `widthMm`. The window grows to the largest stage (unless "size").
 *   A shot may have no `field` (a labelling, arrows or the bare cortex; no `unit` then), and `cortex`, `labelling`, `outlines`, a strip of
 *     several views, `scaleBarKind` "plane", `scaleBarMm` "auto", `scaleBarTile`, `fit`, `dpi`: see panel.ts `Shot`.
 *   A raw shot's `marksInRaster` true draws the bundle's paths, selection edges and points in its pictures too (still reported in the marks);
 *     a view may be a camera direction [x, y, z] (its render `<hemisphere>-from-<x>_<y>_<z>`, `-` as m, `.` as p).
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
    const sh = s.shot, raw = sh.stagePx !== undefined
    if ((sh.field !== undefined && (typeof sh.field !== 'string' || typeof sh.unit !== 'string')) || (!raw && !(sh.widthMm > 0)))
      throw new Error(`${jobPath}: shots[${i}].shot needs field, unit and widthMm (a unit with a field; no field: a labelling, arrows or the bare cortex; widthMm unless stagePx)`)
    if (raw && !(Array.isArray(sh.stagePx) && sh.stagePx.length === 2 && sh.stagePx.every((n) => Number.isInteger(n) && n >= 16 && n <= 8192)))
      throw new Error(`${jobPath}: shots[${i}].shot.stagePx is [width, height], integer pixels 16–8192`)
    if (sh.fit !== undefined && (!raw || !sh.fit || typeof sh.fit !== 'object' || Array.isArray(sh.fit))) throw new Error(`${jobPath}: shots[${i}].shot.fit is {margin?, shared?}, with stagePx`)
    if (sh.marksInRaster !== undefined && (!raw || typeof sh.marksInRaster !== 'boolean')) throw new Error(`${jobPath}: shots[${i}].shot.marksInRaster is true or false, with stagePx`)
    if (sh.views !== undefined && !(Array.isArray(sh.views) && sh.views.every((v) => typeof v === 'string' || (Array.isArray(v) && v.length === 3 && v.every((x) => typeof x === 'number' && Number.isFinite(x))))))
      throw new Error(`${jobPath}: shots[${i}].shot.views is a list of view names or camera directions [x, y, z]`)
  }
  if (j.size !== undefined && !(Array.isArray(j.size) && j.size.length === 2 && j.size.every((n) => Number.isInteger(n) && n > 0))) throw new Error(`${jobPath}: "size" is [width, height] in pixels`)
  if (j.viewer !== undefined) throw new Error(`${jobPath}: "viewer": the viewer is retired; the Editor draws the shots ("editor": its build folder or url)`)
  if (j.editor !== undefined && (typeof j.editor !== 'string' || !j.editor)) throw new Error(`${jobPath}: "editor" is the Editor's build folder or a running Editor's url`)
  const editor = j.editor === undefined ? path.resolve(HERE, '..') : /^https?:/.test(j.editor) ? j.editor : abs(j.editor)
  return { ...j, editor, bundle: abs(j.bundle), out: abs(j.out), tokens: abs(j.tokens), chrome: chromePath(j.chrome) }
}

/** the window: "size" when given, else 1000 × 760, or larger to hold the largest raw stage (`stagePx`) asked */
export function windowOf(job) {
  return job.size ?? job.shots.reduce(([w, h], { shot }) => (shot.stagePx ? [Math.max(w, shot.stagePx[0]), Math.max(h, shot.stagePx[1])] : [w, h]), [1000, 760])
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

/** the release's build (a VERSION file beside index.html: `build: …`), or null for a development build */
export function releaseBuild(dir) {
  const f = path.join(dir, 'VERSION')
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

/** THE EDITOR draws the shots (retire-the-viewer): the page frames the Editor (examples/embed/embed.js), opens the bundle with a scene, sizes the
 *  Stage to each shot's stage and asks it for the panel (`editor.panel`, editor/figure/cortexPanel.ts); its answer is the viewer's `figure` */
async function editorShots(job) {
  if (!fs.existsSync(path.join(job.bundle, 'manifest.json'))) fail(`no bundle at ${job.bundle} (manifest.json)`)
  const dev = /^https?:/.test(job.editor), dist = dev ? null : job.editor
  if (dist && !fs.existsSync(path.join(dist, 'index.html'))) fail(`no Editor build at ${dist} (index.html)`)
  const embedJs = [path.resolve(HERE, '../examples/embed/embed.js'), path.resolve(HERE, 'embed.js')].find((f) => fs.existsSync(f))
  if (!embedJs) fail('embed.js is not beside the driver')
  const tokFile = job.tokens ?? (dist && path.join(dist, 'figure/figure-tokens.json')), tok = tokFile && fs.existsSync(tokFile) ? JSON.parse(fs.readFileSync(tokFile, 'utf8')) : null
  const stageOf = (shot) => shot.stagePx ?? tok?.cortex?.viewer?.size_px ?? [1000, 760]
  const [W0, H0] = windowOf(job), W = Math.max(W0, ...job.shots.map((x) => stageOf(x.shot)[0] + 400)), H = Math.max(H0, ...job.shots.map((x) => stageOf(x.shot)[1] + 400))
  const host = `<!doctype html><meta charset="utf-8"><body style="margin:0"><div id="panel" style="width:${W0}px;height:${H0}px"></div><script type="module">
import { mount } from '/embed.js'
window.__editor = mount(document.getElementById('panel'), { editor: ${JSON.stringify(dev ? job.editor : '/editor/')}, theme: 'light' })
</script>`
  const compute = dist ? [path.join(dist, 'nxr-compute'), path.resolve(HERE, '../../frontend/public/nxr-compute')].find((f) => fs.existsSync(f)) : null
  const figTokens = job.tokens ?? (dist && fs.existsSync(path.join(dist, 'figure/figure-tokens.json')) ? path.join(dist, 'figure/figure-tokens.json') : null)
  server = http.createServer((req, res) => {
    const u = decodeURIComponent(new URL(req.url, 'http://x').pathname)
    if (u === '/' || u === '/__host.html') { res.writeHead(200, { 'content-type': 'text/html' }); res.end(host); return }
    if (u === '/embed.js') { res.writeHead(200, { 'content-type': 'text/javascript' }); res.end(fs.readFileSync(embedJs)); return }
    if (figTokens && u === '/editor/figure/figure-tokens.json') { res.writeHead(200, { 'content-type': 'application/json' }); res.end(fs.readFileSync(figTokens)); return }
    const [root, rel] = u.startsWith('/b/') ? [job.bundle, u.slice(3)] : u.startsWith('/editor/nxr-compute/') && compute ? [compute, u.slice('/editor/nxr-compute/'.length)]
      : u.startsWith('/editor/') && dist ? [dist, u.slice('/editor/'.length) || 'index.html'] : [null, null]
    const f = root && path.resolve(root, rel)
    if (!f || !f.startsWith(path.resolve(root) + path.sep) || !fs.existsSync(f) || !fs.statSync(f).isFile()) { res.writeHead(404); res.end(); return }
    res.writeHead(200, { 'content-type': TYPES[path.extname(f)] ?? 'application/octet-stream', 'access-control-allow-origin': '*' })
    fs.createReadStream(f).pipe(res)
  })
  await new Promise((r) => server.listen(0, '127.0.0.1', r))
  const origin = `http://127.0.0.1:${server.address().port}`
  const { ws, send, evalp, logs } = await devtools(job, W, H)
  await send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: 1, mobile: false })
  await send('Page.navigate', { url: `${origin}/__host.html` })
  const T = 1000 * (job.timeout ?? 120)
  for (let t = Date.now(); !(await evalp('!!window.__editor').catch(() => false)); await sleep(100)) if (Date.now() - t > T) fail('the Editor page did not load\n' + logs.join('\n'))
  const opened = await evalp(`window.__editor.then((e) => e.show({ manifolds: { cortex: { load: { bundle: ${JSON.stringify(origin + '/b/manifest.json')} } } }, layers: [] }))`)
  if (!opened?.ok) fail(`the Editor did not open the bundle: ${opened?.error?.message ?? JSON.stringify(opened)}\n${logs.slice(-20).join('\n')}`)
  await sleep(job.wait ?? 600)
  fs.mkdirSync(job.out, { recursive: true })
  const build = dist ? releaseBuild(dist) : null, done = { ready: { editor: job.editor }, build, shots: [] }
  for (const { name, shot } of job.shots) {
    /* THE STAGE exactly the shot's (the panel holds the Editor's own bars too): measured, the panel resized, until the canvas is W × H */
    const want = stageOf(shot)
    let size = [want[0], want[1]]
    for (let i = 0; i < 6; i++) {
      await evalp(`(() => { const p = document.getElementById('panel'); p.style.width = '${size[0]}px'; p.style.height = '${size[1]}px' })()`)
      await sleep(250)
      const got = (await evalp('window.__editor.then((e) => e.stage())'))?.size
      if (!got) fail('the Editor has no Stage')
      if (got[0] === want[0] && got[1] === want[1]) break
      size = [size[0] + want[0] - got[0], size[1] + want[1] - got[1]]
      if (i === 5) fail(`${name}: the Stage is ${got.join('×')} CSS px, not ${want.join('×')}`)
    }
    let timer
    const head = await Promise.race([evalp(`window.__editor.then(async (e) => { const d = await e.panel(${JSON.stringify(shot)}); if (!d.ok) return { error: d.error?.message ?? 'refused' };`
      + `window.__fig = { svg: d.svg, marks: JSON.stringify(d.marks), report: JSON.stringify(d.report), renders: d.renders || [] };`
      + `return { svg: typeof d.svg === 'string', renders: window.__fig.renders.map((r) => r.label), marked: window.__fig.renders.map((r) => typeof r.marks === 'string'), warnings: d.warnings || [] } })`),
    new Promise((r) => { timer = setTimeout(r, T) }).then(() => fail(`${name}: no figure in ${job.timeout ?? 120} s\n${logs.slice(-20).join('\n')}`))])
    clearTimeout(timer)
    if (!head || head.error) fail(`${name}: the Editor refused the shot: ${head?.error ?? 'no answer'}`)
    await writeShot(job, name, head, evalp, build)
    done.shots.push({ name, renders: head.renders, warnings: head.warnings })
    for (const w of head.warnings) console.error(`snapshot: ${name}: ${w}`)
  }
  fs.writeFileSync(path.join(job.out, 'driver.json'), JSON.stringify(done, null, 2))
  ws.close(); cleanup(); process.exit(0)
}

/** a shot's answer, held by the page in window.__fig, read back in slices and written: <out>/<name>.svg, .marks.json, .report.json, each render */
async function writeShot(job, name, head, evalp, build) {
  const svg = head.svg ? await readSliced(evalp, 'window.__fig.svg') : null, marks = JSON.parse(await readSliced(evalp, 'window.__fig.marks'))
  const report = JSON.parse(await readSliced(evalp, 'window.__fig.report'))
  const pngs = []
  for (const [i, label] of head.renders.entries()) {
    const png = await readSliced(evalp, `window.__fig.renders[${i}].png`)
    if (!png.startsWith('data:image/png;base64,')) fail(`${name}: the render ${label} is no PNG`)
    const own = head.marked[i] ? JSON.parse(await readSliced(evalp, `window.__fig.renders[${i}].marks`)) : null      // raw pictures: each render's marks
    pngs.push([label, Buffer.from(png.slice(png.indexOf(',') + 1), 'base64'), own])
  }
  await evalp('delete window.__fig')
  if (build) report.build = build                                        // a release records its VERSION's build
  const dir = path.join(job.out, name)
  fs.mkdirSync(dir, { recursive: true })
  const sha = job.tokens_sha256 ?? (job.tokens && fs.existsSync(job.tokens) ? createHash('sha256').update(fs.readFileSync(job.tokens)).digest('hex') : undefined)
  const stamp = (m) => (job.tokens ? { ...m, tokens_sha256: sha, tokens_path: job.tokens_path ?? job.tokens } : m)
  for (const [label, buf, own] of pngs) {
    fs.writeFileSync(path.join(dir, `${label}.png`), buf)
    if (own) fs.writeFileSync(path.join(dir, `${label}.marks.json`), JSON.stringify(stamp(own), null, 2))
  }
  if (svg !== null) fs.writeFileSync(path.join(job.out, `${name}.svg`), svg)
  fs.writeFileSync(path.join(job.out, `${name}.marks.json`), JSON.stringify(stamp(marks), null, 2))
  fs.writeFileSync(path.join(job.out, `${name}.report.json`), JSON.stringify(report, null, 2))
}

/** run, not imported: argv[1] is the path as typed, import.meta.url the real one (a library release is a symlink), so both are compared real */
const isMain = () => { try { return fs.realpathSync(process.argv[1]) === fs.realpathSync(fileURLToPath(import.meta.url)) } catch { return false } }
if (process.argv[1] && isMain()) {
  process.on('exit', cleanup)
  const file = process.argv[2]
  if (!file) fail('usage: node snapshot-driver.mjs <job.json>')
  let job
  try { job = parseJob(fs.readFileSync(file, 'utf8'), file) } catch (e) { fail(e.message) }
  if (job.print) await printPanels(job)
  else await editorShots(job)
}
